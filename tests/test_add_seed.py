"""Tests for adding a seed from the page (qaview.project.new_seed, Pipeline.add_seed,
POST /api/seeds): regions read as the page shows them, checks against the FASTA and
the project's seeds, the seed kept in the output folder and read back on open, and
recruited at once with its locus scored next.

Builds on test_rerecruit's synthetic project, with real FASTAs for contig lengths;
recruitments are pickles written by the test, so no recruiter runs.

Run: venv/bin/python tests/test_add_seed.py
"""

import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam

from test_rerecruit import make, put, settle

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def fastas(tmp):
    for name, contig in (("a.fa", "c1"), ("b.fa", "c9")):
        path = os.path.join(tmp, name)
        with open(path, "w") as fh:
            fh.write(f">{contig}\n{'ACGT' * 5000}\n")         # 20 kb
        pysam.faidx(path)


def refused(fn, *args):
    try:
        fn(*args)
    except ValueError as exc:
        return str(exc)
    return None


def test_parse_region():
    from qaview.project import parse_region

    check(parse_region("h1tg000146l:58,730-98,730") == ("h1tg000146l", 58729, 98730),
          "1-based inclusive, commas allowed -> 0-based half-open")
    check(parse_region("c1:1000") == ("c1", 999, 1000), "one position is one base")
    for bad in ("c1", "c1:0-5", "c1:10-5", ""):
        check(refused(parse_region, bad), f"refused: {bad!r}")


def test_checks_save_and_reload():
    from qaview.project import load_project, new_seed, save_added_seed

    with tempfile.TemporaryDirectory() as tmp:
        os.environ["QAVIEW_CACHE"] = tmp   # results beside the test data
        path = make(tmp)
        fastas(tmp)
        p = load_project(path)
        check(refused(new_seed, p, "A", "nope:1-10").startswith("no contig"), "unknown contig")
        check("20,000 bp long" in refused(new_seed, p, "A", "c1:1-20,001"), "past the contig's end")
        check("already a seed: seed1" == refused(new_seed, p, "A", "c1:11-20"),
              "the same region as an existing seed (s.bed's seed1, 10-20 0-based)")
        check("is taken" in refused(new_seed, p, "A", "c1:100-200", "seed2"), "a name in use")
        check("assembly must be" in refused(new_seed, p, "X", "c1:1-10"), "an unknown assembly")
        s = new_seed(p, "B", "c9:1,001-2,000")
        check((s.name, s.assembly, s.contig, s.start, s.end) == ("c9_1000_2000", "B", "c9", 1000, 2000),
              f"a default name, on asm2: {s}")
        save_added_seed(p.output_dir, s)
        again = load_project(path)
        names = [x.name for x in again.seeds]
        check("c9_1000_2000" in names, f"read back when the project is opened: {names}")
        src = [x for x in again.seed_sources if x["format"] == "added in qaview"]
        check(src and src[0]["seeds"] == 1, f"listed as a seed source: {again.seed_sources}")


def test_pipeline_and_endpoint():
    from fastapi.testclient import TestClient

    from qaview.pipeline import Pipeline
    from qaview.project import load_project
    from qaview.server import create_app

    with tempfile.TemporaryDirectory() as tmp:
        os.environ["QAVIEW_CACHE"] = tmp   # results beside the test data
        os.environ["XDG_CONFIG_HOME"] = os.path.join(tmp, "xdg")
        path = make(tmp)
        fastas(tmp)
        p = load_project(path)
        pipe = Pipeline(p)
        put(pipe, "seed1", {}, 100, ["r1"])
        put(pipe, "seed2", {}, 100, ["r1"])
        pipe._recruit_all()
        with TestClient(create_app(p, pipe, start=False)) as c:
            r = c.post("/api/seeds", json={"assembly": "A", "region": "c1:11-20"})
            check(r.status_code == 400 and "already a seed" in r.json()["detail"], f"duplicate: {r.text}")
            r = c.post("/api/seeds", json={"assembly": "A", "region": "c1:x"})
            check(r.status_code == 400 and "contig:start-end" in r.json()["detail"], f"malformed: {r.text}")

            # the recruitment the pipeline will look for, written before the request
            from qaview.project import Seed
            from qatool.SegmentStorage import Segment
            import pickle
            seed = Seed("mine", "added", "A", "c1", 5000, 6000)
            rp = pipe._recruit_path(seed, pipe.params_for("mine"))
            with open(rp, "wb") as fh:
                pickle.dump({"segments": [Segment("A", "c1", 4000, 9000), Segment("B", "c9", 4000, 9000)],
                             "reads": {"r7"}}, fh)
            r = c.post("/api/seeds", json={"assembly": "A", "region": "c1:5,001-6,000", "name": "mine"})
            check(r.status_code == 200 and r.json()["name"] == "mine", f"added: {r.text}")
            check(settle(pipe, ["mine"]), "recruited")
            lid = pipe.seeds["mine"].locus
            check(lid in pipe.loci and pipe.loci[lid].seeds == ["mine"], "its own locus")
            check(pipe._front[:1] == [lid], "scored next")
            check(any(s.name == "mine" for s in p.seeds), "a seed of the project now (stars accept it)")
            r = c.post("/api/stars", json={"seeds": ["mine"], "starred": True})
            check(r.status_code == 200, f"and it can be starred: {r.text}")
            rows = c.get("/api/loci").json()["rows"]
            check(any("mine" in row["seeds"] for row in rows), "its row is in the table")


def main():
    tests = [test_parse_region, test_checks_save_and_reload, test_pipeline_and_endpoint]
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
