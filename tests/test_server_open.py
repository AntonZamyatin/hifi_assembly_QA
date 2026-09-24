"""Tests for opening projects from the page (qaview.server): no project at
start, errors a person can read, and switching projects.

Background work is not started (create_app(start=False)), and the projects are
placeholder files in a temporary folder, so no data is needed. The recent list
goes to a temporary XDG_CONFIG_HOME.

Run: venv/bin/python tests/test_server_open.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fastapi.testclient import TestClient

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


PROJECT = """
name = "{name}"
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
dir = "out_{name}"
"""


def make(tmp, name):
    for f in ("a.fa", "a.bam", "b.fa", "b.bam", "aln.sqlite"):
        open(os.path.join(tmp, f), "w").close()
    with open(os.path.join(tmp, "s.bed"), "w") as fh:
        fh.write("c1\t10\t20\tseed1\n")
    path = os.path.join(tmp, f"{name}.toml")
    with open(path, "w") as fh:
        fh.write(PROJECT.format(name=name))
    return path


def test_open_and_switch():
    from qaview.server import create_app

    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_CONFIG_HOME"] = os.path.join(tmp, "xdg")
        first, second = make(tmp, "first"), make(tmp, "second")
        broken = os.path.join(tmp, "broken.toml")
        with open(broken, "w") as fh:
            fh.write('name = "x"\n[asm1]\nname = "A"\n')

        with TestClient(create_app(start=False)) as c:
            check(c.get("/api/project").json() == {"open": False, "recent": []},
                  "no project at start")
            check(c.get("/api/loci").status_code == 409, "project endpoints wait for a project")

            r = c.post("/api/open", json={"path": os.path.join(tmp, "nope.toml")})
            check(r.status_code == 400 and "no such file" in r.json()["detail"], f"missing: {r.text}")
            r = c.post("/api/open", json={"path": os.path.join(tmp, "s.bed")})
            check(r.status_code == 400 and ".toml" in r.json()["detail"], f"not a toml: {r.text}")
            r = c.post("/api/open", json={"path": broken})
            check(r.status_code == 400 and "missing 'fasta'" in r.json()["detail"],
                  f"a broken project says what is wrong: {r.text}")

            r = c.post("/api/open", json={"path": first})
            check(r.status_code == 200 and r.json()["name"] == "first", f"open: {r.text}")
            p = c.get("/api/project").json()
            check(p["open"] and p["name"] == "first", "the project is open")
            check(c.get("/api/loci").status_code == 200, "and its endpoints answer")

            c.post("/api/open", json={"path": second})
            p = c.get("/api/project").json()
            check(p["name"] == "second", "opening another replaces it")
            check([x["path"] for x in p["recent"]] == [second, first],
                  f"recent, newest first: {p['recent']}")


def main():
    tests = [test_open_and_switch]
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
