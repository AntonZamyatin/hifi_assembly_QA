"""Structural tests for the qaview compute layer (Phase 0).

Synthetic fixtures only — no real BAMs, DBs or FASTAs are read.  The central test is
`test_coverage_matches_reference_*`: the refactored `compute_coverage_bins` must
reproduce the pre-refactor `alignment_plot` binning exactly, so the published figure
is unchanged by the extraction.

Run: venv/bin/python tests/test_compute_layer.py
"""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fixtures import (
    BEST_PID_DELTA, BINSIZE, READS_A, READS_B, SEGS_A, SEGS_B,
    reference_coverage_bins, write_bam, write_db,
)

from qatool.SegmentStorage import Segment
from qatool.plotters.axis import _build_axis_map, plot_to_genome
from qatool.plotters.coverage import alignment_plot, compute_coverage_bins

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def test_coverage_matches_reference_unbucketed(tmp):
    bam = write_bam(os.path.join(tmp, "a.bam"), READS_A)
    cov = compute_coverage_bins(bam, SEGS_A, 0.99, 0.5, BINSIZE, None, BEST_PID_DELTA)
    ref_bins, ref_buckets = reference_coverage_bins(
        bam, SEGS_A, 0.99, 0.5, BINSIZE, None, BEST_PID_DELTA)
    check(cov.bins == ref_bins, f"unbucketed bins differ:\n  new={cov.bins}\n  ref={ref_bins}")
    check(cov.bucket_bins is None and ref_buckets is None, "bucket_bins should be None")


def test_coverage_matches_reference_bucketed(tmp):
    bam = write_bam(os.path.join(tmp, "b.bam"), READS_A)
    buckets = [{"r1", "r2"}, {"r3"}, {"r4"}]     # r5 in no bucket -> falls through to bins
    cov = compute_coverage_bins(bam, SEGS_A, 0.99, 0.5, BINSIZE, buckets, BEST_PID_DELTA)
    ref_bins, ref_buckets = reference_coverage_bins(
        bam, SEGS_A, 0.99, 0.5, BINSIZE, buckets, BEST_PID_DELTA)
    check(cov.bins == ref_bins, f"leftover bins differ:\n  new={cov.bins}\n  ref={ref_bins}")
    check(cov.bucket_bins == ref_buckets, "bucket_bins differ")


def test_alignment_plot_return_unchanged(tmp):
    """The renderer's public return values must survive the extraction."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bam = write_bam(os.path.join(tmp, "c.bam"), READS_A)
    buckets = [{"r1", "r2"}, {"r3"}, {"r4"}]
    fig, ax = plt.subplots()
    _, maxv, bucket_bins = alignment_plot(
        bam, SEGS_A, pid_threshold=0.99, plen_threshold=0.5, binsize=BINSIZE,
        ax=ax, read_buckets=buckets, orientation="x", return_max=True,
        best_pid_delta=BEST_PID_DELTA)
    plt.close(fig)

    ref_bins, ref_buckets = reference_coverage_bins(
        bam, SEGS_A, 0.99, 0.5, BINSIZE, buckets, BEST_PID_DELTA)
    total = np.array(ref_bins, dtype=float)
    for s in ref_buckets:
        total += np.array(s, dtype=float)
    check(bucket_bins == ref_buckets, "alignment_plot bucket_bins changed")
    check(abs(maxv - float(np.max(total))) < 1e-12,
          f"alignment_plot max changed: {maxv} vs {float(np.max(total))}")


def test_bin_genomic_identity(tmp):
    bam = write_bam(os.path.join(tmp, "d.bam"), READS_A)
    cov = compute_coverage_bins(bam, SEGS_A, 0.99, 0.5, BINSIZE, None, BEST_PID_DELTA)
    # 4500bp -> 5 bins, 2200bp -> 3 bins
    check(len(cov.bins) == 8, f"expected 8 bins, got {len(cov.bins)}")
    check(cov.bin_start[:5] == [2000, 3000, 4000, 5000, 6000], f"seg1 starts {cov.bin_start[:5]}")
    check(cov.bin_end[4] == 6500, f"seg1 last bin must clip to 6500, got {cov.bin_end[4]}")
    check(cov.bin_start[5:] == [9000, 10000, 11000], f"seg2 starts {cov.bin_start[5:]}")
    check(cov.bin_end[7] == 11200, f"seg2 last bin must clip to 11200, got {cov.bin_end[7]}")
    check(all(c == "ctgA" for c in cov.bin_contig), "all bins are on ctgA")


def test_bin_read_index(tmp):
    bam = write_bam(os.path.join(tmp, "e.bam"), READS_A)
    cov = compute_coverage_bins(bam, SEGS_A, 0.99, 0.5, BINSIZE, None, BEST_PID_DELTA)
    check(len(cov.bin_offsets) == len(cov.bins) + 1, "CSR offsets must be n_bins+1")
    check(int(cov.bin_offsets[-1]) == len(cov.read_idx), "CSR offsets must end at len(read_idx)")

    # r1 spans 2000-3500 -> bins 0,1 ; r6 was filtered out by pid.
    check("r1" in cov.reads_in_bin(0), f"r1 in bin0, got {cov.reads_in_bin(0)}")
    check("r1" in cov.reads_in_bin(1), f"r1 in bin1, got {cov.reads_in_bin(1)}")
    check("r1" not in cov.reads_in_bin(3), f"r1 not in bin3, got {cov.reads_in_bin(3)}")
    check(all("r6" not in cov.reads_in_bin(i) for i in range(len(cov.bins))),
          "r6 (pid 0.5) must be filtered out everywhere")

    # Every read credited to a bin must have contributed weight to it.
    for i in range(len(cov.bins)):
        if cov.reads_in_bin(i):
            check(cov.bins[i] > 0, f"bin {i} lists reads but has zero coverage")


def test_multimapping_weight(tmp):
    """r5 has two qualifying alignments, so each contributes weight 0.5 per bin.

    Locus 1 (3000-4000) is exactly bin-aligned and covers 1 bin; locus 2
    (9500-10500) straddles and covers 2.  Total = 0.5 * (1 + 2) = 1.5.
    """
    bam = write_bam(os.path.join(tmp, "f.bam"), [r for r in READS_A if r[0] == "r5"])
    cov = compute_coverage_bins(bam, SEGS_A, 0.99, 0.5, BINSIZE, None, BEST_PID_DELTA)
    check(abs(sum(cov.bins) - 1.5) < 1e-12,
          f"expected total weight 1.5 across 3 bins, got {sum(cov.bins)}")
    check(all(v in (0.0, 0.5) for v in cov.bins), f"weights must be 0.5, got {cov.bins}")
    check(sum(1 for v in cov.bins if v > 0) == 3, f"3 bins should be covered, got {cov.bins}")


def test_plot_to_genome_exact():
    segs, total = _build_axis_map(SEGS_A, "asm1")
    check(total == 6700, f"exact-bp total should be 4500+2200=6700, got {total}")
    check(plot_to_genome(segs, 0) == ("ctgA", 2000), "start of segment 1")
    check(plot_to_genome(segs, 4499) == ("ctgA", 6499), "end of segment 1")
    check(plot_to_genome(segs, 4500) == ("ctgA", 9000), "segment 2 starts immediately")
    check(plot_to_genome(segs, 6699) == ("ctgA", 11199), "end of segment 2")
    check(plot_to_genome(segs, 6700) is None, "past the end -> None")
    check(plot_to_genome(segs, -1) is None, "before the start -> None")


def test_plot_to_genome_padded():
    segs, total = _build_axis_map(SEGS_A, "asm1", binsize=BINSIZE)
    check(total == 8000, f"padded total should be 5000+3000=8000, got {total}")
    check(plot_to_genome(segs, 0) == ("ctgA", 2000), "start of segment 1")
    check(plot_to_genome(segs, 4499) == ("ctgA", 6499), "last real base of segment 1")
    check(plot_to_genome(segs, 4600) is None, "inside segment 1's padding -> None")
    check(plot_to_genome(segs, 5000) == ("ctgA", 9000), "segment 2 starts on a bin edge")
    check(plot_to_genome(segs, 7500) is None, "inside segment 2's padding -> None")


def test_region_data_build(tmp):
    from qaview.region import SUPPORT_BUCKET_KEYS, RegionData

    bam1 = write_bam(os.path.join(tmp, "asm1.bam"), READS_A)
    bam2 = write_bam(os.path.join(tmp, "asm2.bam"), READS_B)
    rows = []
    for name, contig, start, length, nm, flag in READS_A:
        rows.append({"assembly": "asm1", "read_name": name, "contig": contig,
                     "start": start, "end": start + length,
                     "aln_type": "secondary" if flag & 256 else "primary"})
    for name, contig, start, length, nm, flag in READS_B:
        rows.append({"assembly": "asm2", "read_name": name, "contig": contig,
                     "start": start, "end": start + length, "aln_type": "primary"})
    db = write_db(os.path.join(tmp, "aln.sqlite"), rows)

    rd = RegionData.build(
        SEGS_A + SEGS_B, {"r1", "r2", "r3"}, db, bam1, bam2,
        asm1="asm1", asm2="asm2", pid_threshold=0.99, plen_threshold=0.5,
        binsize=BINSIZE, best_pid_delta=BEST_PID_DELTA)

    check(rd.plot_h == 6700, f"asm1 extent (exact bp) should be 6700, got {rd.plot_h}")
    check(rd.plot_w == 3300, f"asm2 extent (exact bp) should be 3300, got {rd.plot_w}")
    check(rd.bucket_keys == SUPPORT_BUCKET_KEYS, "bucket key order must mirror report.py")
    check(len(rd.cov1.bins) == 8, f"asm1 has 8 bins, got {len(rd.cov1.bins)}")
    check(len(rd.cov2.bins) == 4, f"asm2 (3300bp) has 4 bins, got {len(rd.cov2.bins)}")

    # Reads aligning to both assemblies (r1, r2, r3) should produce lines.
    check(len(rd.lines) > 0, "expected read-correspondence lines")
    check(len(rd.lines.x0) == len(rd.lines.read_idx) == len(rd.lines.supplementary),
          "line arrays must be parallel")
    check(rd.lines.x0.max() <= rd.plot_w and rd.lines.y0.max() <= rd.plot_h,
          "line endpoints must fall inside the plot extent")
    linked = {rd.read_names[i] for i in rd.lines.read_idx}
    check(linked <= {"r1", "r2", "r3"},
          f"only reads present in both assemblies can be linked, got {linked}")

    check(set(rd.read_class.values()) <= set(SUPPORT_BUCKET_KEYS),
          f"unexpected class labels: {set(rd.read_class.values())}")
    check(rd.genome_at("asm1", 0) == ("ctgA", 2000), "crosshair readout on asm1")
    check(rd.genome_at("asm2", 0) == ("ctgB", 1000), "crosshair readout on asm2")


def test_bin_edges_match_coverage(tmp):
    from qaview.region import RegionData

    bam1 = write_bam(os.path.join(tmp, "g1.bam"), READS_A)
    bam2 = write_bam(os.path.join(tmp, "g2.bam"), READS_B)
    db = write_db(os.path.join(tmp, "g.sqlite"), [
        {"assembly": "asm1", "read_name": "r1", "contig": "ctgA", "start": 2000, "end": 3500},
        {"assembly": "asm2", "read_name": "r1", "contig": "ctgB", "start": 1000, "end": 2500},
    ])
    rd = RegionData.build(SEGS_A + SEGS_B, {"r1"}, db, bam1, bam2,
                          asm1="asm1", asm2="asm2", pid_threshold=0.99,
                          plen_threshold=0.5, binsize=BINSIZE)

    for asm, cov in ((rd.asm1, rd.cov1), (rd.asm2, rd.cov2)):
        left, right = rd.bin_edges(asm)
        check(len(left) == len(cov.bins), f"{asm}: {len(left)} edges vs {len(cov.bins)} bins")
        check(bool(np.all(right > left)), f"{asm}: every bin must have positive width")
        widths = right - left
        check(bool(np.all(widths <= BINSIZE)), f"{asm}: no bin may exceed binsize")
        # Exact-bp edges are contiguous: no padding gaps between segments.
        check(bool(np.all(left[1:] == right[:-1])),
              f"{asm}: exact-bp bin edges must be contiguous, got {left} / {right}")
        check(int(right[-1]) == (rd.plot_h if asm == rd.asm1 else rd.plot_w),
              f"{asm}: last edge must equal the plot extent")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_tests_")
    tests = [
        (test_coverage_matches_reference_unbucketed, True),
        (test_coverage_matches_reference_bucketed, True),
        (test_alignment_plot_return_unchanged, True),
        (test_bin_genomic_identity, True),
        (test_bin_read_index, True),
        (test_multimapping_weight, True),
        (test_plot_to_genome_exact, False),
        (test_plot_to_genome_padded, False),
        (test_region_data_build, True),
        (test_bin_edges_match_coverage, True),
    ]
    for fn, needs_tmp in tests:
        print(f"{fn.__name__} ...")
        fn(tmp) if needs_tmp else fn()

    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
