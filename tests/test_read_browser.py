"""The read browser's queries (region.read_lengths / read_alignments) and the
models the page drives it through.

Run: venv/bin/python tests/test_read_browser.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bokeh.models import Div

from fixtures import make_region_data, write_db

from qaview.glyphs import report_layout
from qaview.region import read_alignments, read_lengths

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _db(tmp):
    rows = [
        # r: two near-best copies and a worse one in A, a supplementary piece, one in B
        dict(assembly="A", read_name="r", contig="c1", start=100, end=10100, pid=0.999, q_len=12000),
        dict(assembly="A", read_name="r", contig="c2", start=500, end=10500, pid=0.9985,
             strand=0, aln_type="secondary", q_len=12000),
        dict(assembly="A", read_name="r", contig="c3", start=0, end=10000, pid=0.995, q_len=12000),
        dict(assembly="A", read_name="r", contig="c1", start=20000, end=22000, pid=0.999,
             plen=0.2, aln_type="supplementary", q_len=12000),
        dict(assembly="B", read_name="r", contig="d1", start=0, end=12000, pid=0.97, q_len=12000),
        dict(assembly="B", read_name="s", contig="d1", start=0, end=8000, q_len=8000),
        dict(assembly="C", read_name="r", contig="x", start=0, end=10, q_len=12000),
    ]
    path = os.path.join(tmp, "rb.sqlite")
    if not os.path.exists(path):
        write_db(path, rows)
    return path


def test_lengths(tmp):
    db = _db(tmp)
    got = read_lengths(db, ["r", "s", "nope"], ["A", "B"])
    check(got == {"r": 12000, "s": 8000}, f"lengths {got}")


def test_alignments_marked(tmp):
    db = _db(tmp)
    rows = read_alignments(db, "r", ["A", "B"], pid=0.997, plen=0.99, tol=0.001)
    key = {(a["assembly"], a["contig"], a["start"]): a for a in rows}
    check(len(rows) == 5, f"every alignment in A and B, none in C: {len(rows)}")
    check(key[("A", "c1", 100)]["best"] and key[("A", "c2", 500)]["best"],
          "both copies within best_pid_delta are best placements")
    check(key[("A", "c2", 500)]["strand"] == "-", "strand follows the record")
    c3 = key[("A", "c3", 0)]
    check(not c3["passes"] and not c3["best"], "below the display pid: listed, unmarked")
    sup = key[("A", "c1", 20000)]
    check(sup["type"] == "supplementary" and not sup["passes"], "a short piece fails plen")
    check(not key[("B", "d1", 0)]["passes"], "B's alignment is below the display pid")


def test_page_models(tmp):
    rd = make_region_data(tmp, "pm")
    layout = report_layout(rd)
    picks = [c for c in layout.children if isinstance(c, Div) and c.name == "qv-pick"]
    check(len(picks) == 1 and not picks[0].visible, "a hidden qv-pick in the layout")
    check(len(picks[0].js_property_callbacks.get("change:text", [])) == 1,
          "setting its text selects the read")
    grid = [c for c in layout.children if not isinstance(c, Div)][0]
    taps = [cb for child, _r, _c in grid.children if child is not None
            for cb in child.js_event_callbacks.get("tap", [])]
    check(len(taps) == 2 and all("qv-bin" in cb.code for cb in taps),
          "a click on either marginal sends its bin's reads")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_rb_tests_")
    tests = [test_lengths, test_alignments_marked, test_page_models]
    for fn in tests:
        print(f"{fn.__name__} ...")
        fn(tmp)
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
