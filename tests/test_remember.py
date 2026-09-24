"""Tests for remembered alignments (qaview.remember), in a temporary folder.

Run: venv/bin/python tests/test_remember.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qaview import remember

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def alignment(strand="+"):
    return {"top": "  ACGT-ACGT", "mid": "   ||| ||||", "bottom": "GGACGTTACGT", "strand": strand,
            "top_name": "A", "top_contig": "a1", "top_start": 1000, "top_end": 1008,
            "bottom_name": "B", "bottom_contig": "b1", "bottom_start": 5000, "bottom_end": 5011,
            "compared": 9, "identity": 1.0, "mismatches": 0, "ins": 0, "ins_opens": 0,
            "del": 1, "del_opens": 1, "outside_top": 0, "outside_bottom": 2}


def test_save_load_list_forget():
    with tempfile.TemporaryDirectory() as out:
        m = remember.save(out, alignment(), note="look here", seeds=["s/1", "s2"], locus="L",
                          settings={"pid": 0.997})
        check(os.path.exists(os.path.join(out, "remembered", m["id"] + ".json")), "json written")
        check(os.path.exists(os.path.join(out, "remembered", m["id"] + ".txt")), "text written")
        check("/" not in m["id"] and m["id"].startswith("s_1_"), f"id is a safe file name: {m['id']}")
        e = remember.load(out, m["id"])
        check(e["alignment"]["top"] == alignment()["top"] and e["note"] == "look here",
              "the snapshot comes back as saved")
        check([x["id"] for x in remember.all_meta(out)] == [m["id"]] and
              "alignment" not in remember.all_meta(out)[0], "listing carries no alignment")
        m2 = remember.save(out, alignment(), note="", seeds=["s/1"], locus="L", settings={})
        check(m2["id"] != m["id"], "two in the same second get distinct ids")
        remember.forget(out, m["id"])
        check([x["id"] for x in remember.all_meta(out)] == [m2["id"]], "forget removes one")
        try:
            remember.save(out, {"top": "A"}, note="", seeds=["s"], locus="L", settings={})
            check(False, "an incomplete alignment is refused")
        except ValueError:
            pass


def test_place_by_genome_coordinates():
    meta = {"id": "x", "note": "", "created": "t",
            "top": {"assembly": "A", "contig": "a1", "start": 1000, "end": 1008},
            "bottom": {"assembly": "B", "contig": "b1", "start": 5000, "end": 5011, "strand": "+"}}
    s1 = [("a0", 0, 500, 0), ("a1", 900, 2000, 500)]      # (contig, start, end, offset)
    s2 = [("b1", 4000, 6000, 0)]
    got = remember.place([meta], s1, s2, "A", "B")
    check(got == [{"id": "x", "note": "", "created": "t", "x0": 1000, "x1": 1011,
                   "y0": 600, "y1": 608}], f"rectangle in plot coordinates: {got}")
    check(remember.place([meta], s1, [("b2", 0, 9000, 0)], "A", "B") == [],
          "not drawn where the stretch is not in the locus")


def test_text_alignment_coordinates():
    for strand, first_bottom in (("+", "5,001"), ("-", "5,011")):
        e = {"id": "i", "created": "t", "note": "", "seeds": ["s"], "settings": {},
             "top": {"assembly": "A", "contig": "a1", "start": 1000, "end": 1008},
             "bottom": {"assembly": "B", "contig": "b1", "start": 5000, "end": 5011, "strand": strand},
             "alignment": alignment(strand)}
        rows = [l for l in remember.text_alignment(e).splitlines() if l and not l.startswith("#")]
        check(rows[0].split()[:2] == ["A", "1,001"] and rows[2].split()[:2] == ["B", first_bottom],
              f"strand {strand}: first bases at {rows[0].split()[1]} and {rows[2].split()[1]}")


def main():
    tests = [test_save_load_list_forget, test_place_by_genome_coordinates,
             test_text_alignment_coordinates]
    for fn in tests:
        print(f"{fn.__name__} ...")
        fn()
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
