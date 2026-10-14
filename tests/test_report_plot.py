"""report_plot, the notebooks' figure, end to end on fixtures (verbose, so the
likelihood report and its DB queries run too).

Run: venv/bin/python tests/test_report_plot.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import pysam

from fixtures import READS_A, READS_B, SEGS_A, SEGS_B, write_bam, write_db

from qatool.SegmentStorage import SegmentCollection
from qatool.plotters import report_plot

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def test_renders(tmp):
    bam1 = write_bam(os.path.join(tmp, "e1.bam"), READS_A)
    bam2 = write_bam(os.path.join(tmp, "e2.bam"), READS_B)
    db = write_db(os.path.join(tmp, "e.sqlite"), [
        {"assembly": "asm1", "read_name": "r1", "contig": "ctgA", "start": 2000, "end": 3500},
        {"assembly": "asm2", "read_name": "r1", "contig": "ctgB", "start": 1000, "end": 2500},
    ])
    fa = os.path.join(tmp, "e.fa")
    with open(fa, "w") as fh:
        fh.write(">ctgA\n" + "ACGTTGCA" * 2500 + "\n>ctgB\n" + "ACGTTGCA" * 1500 + "\n")
    pysam.faidx(fa)
    out = os.path.join(tmp, "fig.png")
    report_plot("", "", SegmentCollection([SEGS_A[0], SEGS_B[0]]), {"r1"}, db, 0.99, 0.5,
                output_path=out, dpi=60, asm1_bam=bam1, asm2_bam=bam2,
                asm1_fasta=fa, asm2_fasta=fa, asm1="asm1", asm2="asm2",
                pid_tie_eps=0.0001, best_pid_delta=0.001, show=False, verbose=True)
    check(os.path.exists(out) and os.path.getsize(out) > 5_000, "a real figure is written")


def main():
    tmp = tempfile.mkdtemp(prefix="report_plot_tests_")
    tests = [test_renders]
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
