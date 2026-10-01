"""Tests for recruiting seeds again with parameters of their own (qaview.pipeline,
the /api/recruit endpoints): seeds move between loci, a locus with no seed left
goes, parameters are kept per seed across a reopen, and failures leave the seed
where it was.

Recruitments are pickles written by the test under the keys the pipeline asks
for, so no data or recruiter run is needed; the failure case runs one real
worker on placeholder files. Scoring is never started.

Run: venv/bin/python tests/test_rerecruit.py
"""

import os
import pickle
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


PROJECT = """
name = "t"
[asm1]
name = "A"
fasta = "a.fa"
bam = "a.bam"
[asm2]
name = "B"
fasta = "b.fa"
bam = "b.bam"
[data]
db = "aln.sqlite"
[[seeds]]
format = "bed"
path = "s.bed"
[output]
dir = "out"
"""


def make(tmp):
    for f in ("a.fa", "a.bam", "b.fa", "b.bam", "aln.sqlite"):
        open(os.path.join(tmp, f), "w").close()
    with open(os.path.join(tmp, "s.bed"), "w") as fh:
        fh.write("c1\t10\t20\tseed1\nc1\t30\t40\tseed2\n")
    path = os.path.join(tmp, "t.toml")
    with open(path, "w") as fh:
        fh.write(PROJECT)
    return path


def put(pipe, name, params, end, reads):
    """Write the recruitment `pipe` would make for this seed and these parameters."""
    from qatool.SegmentStorage import Segment

    st = pipe.seeds[name]
    path = pipe._recruit_path(st.seed, {**pipe.project.recruitment, **params})
    with open(path, "wb") as fh:
        pickle.dump({"segments": [Segment("A", "c1", 0, end), Segment("B", "c9", 0, end)],
                     "reads": set(reads)}, fh)


def settle(pipe, names, timeout=60.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        with pipe._lock:
            if all(pipe.seeds[n].pending is None for n in names):
                return True
        time.sleep(0.05)
    return False


def locus_of(pipe, name):
    return pipe.seeds[name].locus


def test_rerecruit():
    from qaview.pipeline import Pipeline
    from qaview.project import load_project

    with tempfile.TemporaryDirectory() as tmp:
        project = load_project(make(tmp))
        pipe = Pipeline(project)
        defaults = dict(project.recruitment)
        put(pipe, "seed1", {}, 100, ["r1", "r2"])
        put(pipe, "seed2", {}, 100, ["r1", "r2"])
        pipe._recruit_all()                               # both on disk: no workers
        check(locus_of(pipe, "seed1") == locus_of(pipe, "seed2") and len(pipe.loci) == 1,
              "two seeds recruiting the same locus share it")
        first = locus_of(pipe, "seed1")
        seq0, _ = pipe.changes_since(-1)

        # seed1 with a smaller merge_gap recruits a smaller locus
        own = {"merge_gap": 1000}
        put(pipe, "seed1", own, 50, ["r1"])
        out = pipe.rerecruit(["seed1"], {**defaults, **own})
        check(out["queued"] == ["seed1"] and out["unchanged"] == [], f"queued: {out}")
        check(settle(pipe, ["seed1"]), "the re-recruitment finishes")
        second = locus_of(pipe, "seed1")
        check(second != first and second in pipe.loci, "seed1 moved to a new locus")
        check(pipe.loci[first].seeds == ["seed2"], "seed2 stays in the old one")
        check(pipe.loci[first].recruiting == set(), "the old locus is no longer busy")
        check(pipe.overrides.get("seed1") == own, f"only what differs is kept: {pipe.overrides.get('seed1')}")
        check(pipe.overrides.get("seed2") == {}, "seed2 has none")
        check(pipe._front[:1] == [second], "the new locus is scored next")
        rows = {r["id"]: r for r in pipe.rows()}
        check(rows[second]["custom"] == ["seed1"] and rows[first]["custom"] == [],
              "rows say which seeds have parameters of their own")
        info = pipe.seed_info(second)
        check(info[0]["params"]["merge_gap"] == 1000 and info[0]["own"] == own,
              f"seed_info: {info}")

        out = pipe.rerecruit(["seed1"], {**defaults, **own})
        check(out["queued"] == [] and out["unchanged"] == ["seed1"], "the same parameters again: nothing to do")

        # seed2 follows: the old locus has no seed left and goes
        put(pipe, "seed2", own, 50, ["r1"])
        pipe.rerecruit(["seed2"], {**defaults, **own})
        check(settle(pipe, ["seed2"]), "seed2 finishes")
        check(first not in pipe.loci, "a locus with no seed left is removed")
        check(pipe.loci[second].seeds == ["seed1", "seed2"], "both seeds share the new locus")
        _, changed = pipe.changes_since(seq0)
        gone = [r for r in changed if r["id"] == first]
        check(gone == [{"id": first, "removed": True, "moved_to": second}],
              f"the removal is reported with where the seeds went: {gone}")

        # the parameters are kept for a reopened project
        again = Pipeline(project)
        check(again.seeds["seed1"].params["merge_gap"] == 1000
              and again.seeds["seed1"].path == pipe.seeds["seed1"].path,
              "a reopened project recruits seed1 with its own parameters")
        again._recruit_all()
        check(set(again.loci) == {second}, "and finds the same loci on disk")

        # back to the project file: the old recruitment is on disk and the old locus returns
        pipe.rerecruit(["seed1"], defaults)
        check(settle(pipe, ["seed1"]), "reverting finishes")
        check(locus_of(pipe, "seed1") == first and first in pipe.loci, "the first locus is back")
        check(pipe.overrides.get("seed1") == {}, "and seed1's own parameters are cleared")
        _, changed = pipe.changes_since(seq0)
        check(any(r["id"] == first and not r.get("removed") for r in changed),
              "the returning locus is reported as a row again")

        # bad requests
        for bad, what in (({"merge_gap": "x"}, "whole number"), ({"pid_threshold": 2}, "between"),
                          ({"nope": 1}, "unknown")):
            try:
                pipe.rerecruit(["seed1"], {**defaults, **bad})
                check(False, f"{bad} is refused")
            except ValueError as exc:
                check(what in str(exc), f"{bad}: {exc}")
        try:
            pipe.rerecruit(["nobody"], defaults)
            check(False, "an unknown seed is refused")
        except KeyError:
            pass


def test_failure_keeps_the_seed():
    from qaview.pipeline import Pipeline
    from qaview.project import load_project

    with tempfile.TemporaryDirectory() as tmp:
        project = load_project(make(tmp))
        pipe = Pipeline(project)
        put(pipe, "seed1", {}, 100, ["r1"])
        put(pipe, "seed2", {}, 100, ["r1"])
        pipe._recruit_all()
        lid = locus_of(pipe, "seed1")
        try:
            # nothing on disk for these parameters: a real worker runs, on empty files
            pipe.rerecruit(["seed1"], {**project.recruitment, "ext_rounds": 7})
            check(settle(pipe, ["seed1"], timeout=120), "the failing worker returns")
        finally:
            pipe.stop()
        loc = pipe.loci[lid]
        check(locus_of(pipe, "seed1") == lid and "seed1" in loc.seeds, "seed1 keeps its locus")
        check(loc.recruiting == set() and (loc.recruit_error or "").startswith("seed1:"),
              f"the row carries the error: {loc.recruit_error}")
        check(pipe.overrides.get("seed1") == {}, "a failed request does not keep its parameters")
        check(pipe.seeds["seed1"].params == project.recruitment, "nor change the seed's")


def test_endpoints():
    from fastapi.testclient import TestClient

    from qaview.pipeline import Pipeline
    from qaview.project import load_project
    from qaview.server import create_app

    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_CONFIG_HOME"] = os.path.join(tmp, "xdg")
        project = load_project(make(tmp))
        pipe = Pipeline(project)
        put(pipe, "seed1", {}, 100, ["r1"])
        put(pipe, "seed2", {}, 100, ["r1"])
        pipe._recruit_all()
        lid = locus_of(pipe, "seed1")
        with TestClient(create_app(project, pipe, start=False)) as c:
            ranges = c.get("/api/project").json()["recruitment_ranges"]
            check(ranges["merge_gap"][0] == "int" and ranges["pid_threshold"][0] == "float",
                  f"the page learns each parameter's type: {ranges}")
            d = c.get(f"/api/loci/{lid}/seeds").json()
            check([s["name"] for s in d["seeds"]] == ["seed1", "seed2"]
                  and d["seeds"][0]["start"] == 10 and d["defaults"] == project.recruitment,
                  f"seeds: {d}")
            check(c.get("/api/loci/nope/seeds").status_code == 404, "unknown locus")
            params = {k: str(v) for k, v in project.recruitment.items()}   # as the form sends them
            r = c.post("/api/recruit", json={"seeds": ["seed1"], "params": params})
            check(r.status_code == 200 and r.json()["unchanged"] == ["seed1"] and "seq" in r.json(),
                  f"strings from the form are read as numbers: {r.text}")
            r = c.post("/api/recruit", json={"seeds": ["seed1"], "params": {**params, "merge_gap": "1.5"}})
            check(r.status_code == 400 and "whole number" in r.json()["detail"], f"bad value: {r.text}")
            r = c.post("/api/recruit", json={"seeds": ["x"], "params": params})
            check(r.status_code == 400 and "unknown seeds" in r.json()["detail"], f"unknown: {r.text}")
            r = c.post("/api/recruit", json={"seeds": [], "params": params})
            check(r.status_code == 400, "no seeds")


def main():
    tests = [test_rerecruit, test_failure_keeps_the_seed, test_endpoints]
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
