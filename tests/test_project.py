"""Tests for the qaview project file and seed sources (docs/qaview.md §13.2).

Builds small projects in a temporary directory, so no real data is needed: the
files named are empty placeholders, which is all the loader checks for.

Run: venv/bin/python tests/test_project.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qaview.project import (RECRUITMENT_DEFAULTS, ProjectError, load_bed, load_project,
                            load_quast)

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


QUAST = ("contig\tstart\tend\ttype\tinconsistency_bp\tbp_start\tbp_end\n"
         "c1\t100\t500\trelocation\t\t300\t301\n"
         "c1\t900\t1500\trelocation\t12\t1200\t1250\n"
         "c2\t0\t400\tinversion\t\t200\t201\n")


def make(tmp, toml_body, files=("a.fa", "a.bam", "b.fa", "b.bam", "aln.sqlite")):
    for f in files:
        open(os.path.join(tmp, f), "w").close()
    with open(os.path.join(tmp, "q.tsv"), "w") as fh:
        fh.write(QUAST)
    path = os.path.join(tmp, "p.toml")
    with open(path, "w") as fh:
        fh.write(toml_body)
    return path


BASE = """
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
[output]
dir = "out"
"""


def test_quast_names_match_the_notebooks():
    with tempfile.TemporaryDirectory() as tmp:
        make(tmp, "")
        seeds = load_quast(os.path.join(tmp, "q.tsv"), "A")
    check([s.name for s in seeds] == ["c1_relocation_001", "c1_relocation_002",
                                      "c2_inversion_001"], f"names: {[s.name for s in seeds]}")
    check(seeds[1].meta == {"type": "relocation", "breakpoint": "1,200-1,250",
                            "inconsistency_bp": 12}, f"meta: {seeds[1].meta}")
    check(seeds[0].interval == ("A", "c1", 100, 500), "interval")


def test_bed_names_and_skips():
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "s.bed")
        with open(p, "w") as fh:
            fh.write("track name=x\n# comment\n\nc1\t10\t20\tfirst\nc2\t5\t9\n")
        seeds = load_bed(p, "B")
        check([s.name for s in seeds] == ["first", "c2_5_9"], f"names: {[s.name for s in seeds]}")
        with open(p, "w") as fh:
            fh.write("c1\t20\t20\n")
        try:
            load_bed(p, "B")
            check(False, "an empty interval must be rejected")
        except ProjectError:
            pass


def test_project_defaults_and_sources():
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "s.bed"), "w") as fh:
            fh.write("c1\t100\t500\tsame_as_quast\nc9\t1\t50\n")
        path = make(tmp, BASE + '[[seeds]]\nformat = "quast"\npath = "q.tsv"\n'
                              '[[seeds]]\nformat = "bed"\npath = "s.bed"\nassembly = "A"\n')
        p = load_project(path)
        check(p.library_assembly == "B", "lambda defaults to asm2")
        check(p.recruitment == RECRUITMENT_DEFAULTS, "recruitment defaults to the notebooks'")
        check(p.resources == {"recruit_workers": 32, "score_workers": 32}, "resources default 32")
        check(p.output_dir == os.path.join(tmp, "out"), "paths resolve against the file")
        check(len(p.seeds) == 4, f"3 QUAST + 2 BED with one shared interval -> 4, got {len(p.seeds)}")
        first = p.seeds[0]
        check(first.name == "c1_relocation_001" and first.meta.get("aliases") == ["same_as_quast"],
              f"a shared interval is one seed with an alias: {first}")
        check([s["seeds"] for s in p.seed_sources] == [3, 2], "per-source counts")
        os.environ["QAVIEW_OUTPUT"] = os.path.join(tmp, "mine")
        try:
            p = load_project(path)
        finally:
            del os.environ["QAVIEW_OUTPUT"]
        check(p.output_dir == os.path.join(tmp, "mine", p.name), "--output: ROOT/<project name>")


def test_errors_are_reported():
    with tempfile.TemporaryDirectory() as tmp:
        path = make(tmp, BASE + '[[seeds]]\nformat = "quast"\npath = "q.tsv"\n',
                    files=("a.fa", "b.fa", "b.bam", "aln.sqlite"))
        try:
            load_project(path)
            check(False, "a missing BAM must be reported")
        except ProjectError as exc:
            check("a.bam" in str(exc), f"message names the file: {exc}")
        for body, what in [
            (BASE + '[[seeds]]\nformat = "vcf"\npath = "q.tsv"\n', "unknown seed format"),
            (BASE + '[[seeds]]\nformat = "quast"\npath = "q.tsv"\nassembly = "C"\n',
             "seed assembly that is neither candidate"),
            (BASE + '[recruitment]\npid = 0.99\n[[seeds]]\nformat = "quast"\npath = "q.tsv"\n',
             "unknown recruitment key"),
            (BASE + '[recruitment]\npid_threshold = "high"\n[[seeds]]\nformat = "quast"\npath = "q.tsv"\n',
             "a recruitment value that is not a number"),
            (BASE + '[recruitment]\nmerge_gap = -1\n[[seeds]]\nformat = "quast"\npath = "q.tsv"\n',
             "a recruitment value out of range"),
            (BASE, "no seed sources"),
        ]:
            path = make(tmp, body)
            try:
                load_project(path)
                check(False, f"{what} must be rejected")
            except ProjectError:
                pass


def test_example_project_parses():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = load_project(os.path.join(root, "examples/projects/hg002_hic.toml"), check_files=False)
    check((p.asm1.name, p.asm2.name) == ("hifiasm_hic", "curated"), "example candidates")
    check(p.loci_dir.endswith("examples/analysis/hg002_hic/loci"),
          "the example shares the notebooks' score cache")


def main():
    tests = [test_quast_names_match_the_notebooks, test_bed_names_and_skips,
             test_project_defaults_and_sources, test_errors_are_reported,
             test_example_project_parses]
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
