"""Seed picking from a QUAST misassembly table (qaview.seeds + the sidebar table).

Synthetic fixtures only.  Pins the region naming shared with the analysis
notebooks, the optional notebook score, and that picking a row fills the seed.

Run: venv/bin/python tests/test_seeds.py
"""

import json
import os
import pickle
import sys
import tempfile

from qaview.app import QaviewApp
from qaview.recruit import cache_path
from qaview.seeds import load_seeds, scores_dir, seeds_path

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


TSV = """contig\tstart\tend\ttype\tinconsistency_bp\tbp_start\tbp_end
c1\t1000\t41000\ttranslocation\t\t21000\t21001
c1\t5000\t45000\ttranslocation\t\t25000\t25400
c1\t9000\t49000\trelocation\t684\t29000\t29001
c2\t100\t40100\ttranslocation\t\t20100\t20101
"""


def make_base(tmp, run="r"):
    base = os.path.join(tmp, "base")
    os.makedirs(os.path.dirname(seeds_path(base, run)))
    with open(seeds_path(base, run), "w") as fh:
        fh.write(TSV)
    os.makedirs(scores_dir(base, run))
    with open(os.path.join(scores_dir(base, run), "c1_translocation_002.json"), "w") as fh:
        json.dump({"delta_s": -12.345}, fh)
    return base


def test_region_names_match_notebooks(tmp):
    base = make_base(tmp)
    df = load_seeds(seeds_path(base, "r"), scores_dir(base, "r"))
    check(list(df["region"]) == ["c1_translocation_001", "c1_translocation_002",
                                 "c1_relocation_001", "c2_translocation_001"],
          f"numbering counts rows per (contig, type) in file order: {list(df['region'])}")
    check(df["breakpoint"][0] == "21,000" and df["breakpoint"][1] == "25,000-25,400",
          f"a 1 bp breakpoint is a point, a longer one a range: {list(df['breakpoint'])}")
    check(df["inconsistency_bp"][2] == 684, "inconsistency carried through")
    check(df["delta_s"][1] == -12.35 and df["delta_s"].isna().sum() == 3,
          "the notebook's Delta S is attached where its JSON exists, missing elsewhere")


def test_pick_fills_seed(tmp):
    base = make_base(tmp)
    app = QaviewApp(cache_dir=os.path.join(tmp, "cache"), base=base, quast_run="r")
    check(len(app.seed_table.value) == 4, "table loads from <base>/quast/<run>")
    calls = []
    app._on_recruit = lambda *_: calls.append(1)
    app.seed_table.selection = [2]
    check((app.contig.value, app.start.value, app.end.value) == ("c1", 9000, 49000),
          "picking a row sets contig/start/end")
    check(not calls, "an uncached seed must not start a minutes-long recruitment")


def test_cached_seed_recruits_at_once(tmp):
    base = make_base(tmp)
    app = QaviewApp(cache_dir=os.path.join(tmp, "cache"), base=base, quast_run="r")
    params = {**app._params(), "contig": "c2", "start": 100, "end": 40100}
    os.makedirs(app.cache_dir, exist_ok=True)
    with open(cache_path(app.cache_dir, params), "wb") as fh:
        pickle.dump({}, fh)
    app._load_seeds()
    check(list(app.seed_table.value["cached"]) == ["", "", "", "✓"],
          f"cached flag follows the recruitment cache: {list(app.seed_table.value['cached'])}")
    calls = []
    app._on_recruit = lambda *_: calls.append(1)
    app.seed_table.selection = [3]
    check(calls == [1], "a cached seed loads immediately")


def test_missing_table(tmp):
    app = QaviewApp(cache_dir=os.path.join(tmp, "cache"), base=os.path.join(tmp, "none"))
    check(len(app.seed_table.value) == 0 and "no seed table" in app.seed_note.object,
          "a missing table empties the list and says where it looked")


def main():
    tests = [test_region_names_match_notebooks, test_pick_fills_seed,
             test_cached_seed_recruits_at_once, test_missing_table]
    for fn in tests:
        print(f"{fn.__name__} ...")
        with tempfile.TemporaryDirectory() as tmp:
            fn(tmp)
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
