"""Tests for boundary clipping of read-alignment blocks.

Blocks that overhang a displayed segment used to be discarded whole, losing
2L/(S+L) of reads — ~49% on real 41.5 kb segments, 100% once a segment is shorter
than a read.  They are now clipped to the visible part.  The property that matters is
that clipping happens in *read* space: clipping each axis independently would leave a
line's two ends describing different stretches of read, silently faking the slope.

Run: venv/bin/python tests/test_clipping.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

from fixtures import write_bam, write_db

from qatool.SegmentStorage import Segment
from qatool.plotters.reads_dotplot import (
    _map_in_segment, _overlapping_segment, _read_to_ref, _ref_to_read, _visible_read_span,
)

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


LOOKUP = {"ctgA": [(2000, 6500, 0)]}          # one displayed segment, plot offset 0


def test_ref_to_read_inverts():
    for strand in (0, 1):
        for rp in (0, 250, 1500, 2999):
            ref = _read_to_ref(rp, 0, 3000, 5000, 8000, strand)
            back = _ref_to_read(ref, 0, 3000, 5000, 8000, strand)
            check(abs(back - rp) < 1e-6,
                  f"strand {strand}: read {rp} -> ref {ref} -> read {back}")


def test_overlapping_segment_picks_largest():
    lk = {"ctgA": [(0, 1000, 0), (1000, 9000, 1000)]}
    check(_overlapping_segment(lk, "ctgA", 900, 5000) == (1000, 9000, 1000),
          "must pick the segment with the greatest overlap, not the first")
    check(_overlapping_segment(lk, "ctgA", 20_000, 30_000) is None, "no overlap -> None")
    check(_overlapping_segment(lk, "ctgB", 0, 100) is None, "unknown contig -> None")


def test_visible_span_forward():
    """Block ref 5000-8000 overhangs the segment end (6500) by 1500 of its 3000 bp."""
    v = _visible_read_span(LOOKUP, "ctgA", 5000, 8000, 0, 3000, strand=1)
    check(v is not None, "overhanging block must now be visible, not dropped")
    lo, hi, seg = v
    check(abs(lo - 0) < 1e-6 and abs(hi - 1500) < 1e-6,
          f"visible read span should be 0..1500, got {lo}..{hi}")
    check(seg == (2000, 6500, 0), f"wrong segment {seg}")


def test_visible_span_reverse():
    """On the reverse strand the read's 5' end maps to high reference coordinates,
    so clipping the high-ref overhang removes the START of the read span."""
    v = _visible_read_span(LOOKUP, "ctgA", 5000, 8000, 0, 3000, strand=0)
    lo, hi, _ = v
    check(abs(lo - 1500) < 1e-6 and abs(hi - 3000) < 1e-6,
          f"reverse-strand visible span should be 1500..3000, got {lo}..{hi}")


def test_fully_contained_unchanged():
    v = _visible_read_span(LOOKUP, "ctgA", 3000, 4000, 0, 1000, strand=1)
    lo, hi, _ = v
    check(abs(lo - 0) < 1e-6 and abs(hi - 1000) < 1e-6,
          f"a contained block must be untouched, got {lo}..{hi}")


def test_block_outside_still_dropped():
    check(_visible_read_span(LOOKUP, "ctgA", 8000, 9000, 0, 1000, 1) is None,
          "a block with no overlap must still be dropped, not clamped to an edge")


def test_map_in_segment_clamps():
    seg = (2000, 6500, 0)
    check(_map_in_segment(seg, 2000) == 0, "segment start maps to its offset")
    check(_map_in_segment(seg, 6500) == 4500, "segment end maps to offset + length")
    check(_map_in_segment(seg, 9999) == 4500, "float overshoot must clamp, not run off-axis")
    check(_map_in_segment(seg, 0) == 0, "undershoot must clamp too")


def test_clipped_line_keeps_slope(tmp):
    """The clipped line must lie on the same diagonal the unclipped one would."""
    from qaview.region import RegionData

    segs = [Segment("asm1", "ctgA", 2_000, 6_500), Segment("asm2", "ctgB", 1_000, 4_300)]
    # 'ov' overhangs BOTH segments; under the old rule it produced no line at all.
    reads_a = [("ov", "ctgA", 5_000, 3_000, 0, 0)]
    reads_b = [("ov", "ctgB", 2_000, 3_000, 0, 0)]
    bam1 = write_bam(os.path.join(tmp, "c1.bam"), reads_a)
    bam2 = write_bam(os.path.join(tmp, "c2.bam"), reads_b)
    db = write_db(os.path.join(tmp, "c.sqlite"), [
        {"assembly": "asm1", "read_name": "ov", "contig": "ctgA", "start": 5_000, "end": 8_000},
        {"assembly": "asm2", "read_name": "ov", "contig": "ctgB", "start": 2_000, "end": 5_000},
    ])

    rd = RegionData.build(segs, {"ov"}, db, bam1, bam2, asm1="asm1", asm2="asm2",
                          pid_threshold=0.99, plen_threshold=0.5, binsize=1000)
    check(len(rd.lines) == 1, f"expected one clipped line, got {len(rd.lines)}")
    if not len(rd.lines):
        return
    x0, y0 = rd.lines.x0[0], rd.lines.y0[0]
    x1, y1 = rd.lines.x1[0], rd.lines.y1[0]

    check(0 <= min(x0, x1) and max(x0, x1) <= rd.plot_w, f"x {x0},{x1} outside 0..{rd.plot_w}")
    check(0 <= min(y0, y1) and max(y0, y1) <= rd.plot_h, f"y {y0},{y1} outside 0..{rd.plot_h}")

    # Both alignments advance with the read at 1 ref bp per read bp, so the visible
    # correspondence must sit exactly on a 45-degree diagonal.
    dx, dy = x1 - x0, y1 - y0
    check(abs(dx - dy) < 1e-6, f"slope distorted by clipping: dx={dx}, dy={dy}")

    # asm1 is clipped at ref 6500 (1500 of 3000 bp visible); asm2 at ref 4300
    # (2300 visible).  The shared visible read span is therefore 1500 bp.
    check(abs(abs(dx) - 1500) < 1e-6,
          f"visible span should be the 1500 bp both sides share, got {abs(dx)}")


def test_overlay_recovers_overhanging_read(tmp):
    """report_plot's overlay path must recover the same read."""
    from matplotlib import pyplot as plt
    from matplotlib.collections import LineCollection
    from qatool.plotters.reads_overlay import reads_alignment_dotplot

    segs = [Segment("asm1", "ctgA", 2_000, 6_500), Segment("asm2", "ctgB", 1_000, 4_300)]
    db = write_db(os.path.join(tmp, "o.sqlite"), [
        {"assembly": "asm1", "read_name": "ov", "contig": "ctgA",
         "start": 5_000, "end": 8_000, "pid": 0.999},
        {"assembly": "asm2", "read_name": "ov", "contig": "ctgB",
         "start": 2_000, "end": 5_000, "pid": 1.0},
    ])
    fig, ax = plt.subplots()
    reads_alignment_dotplot(segs, {"ov"}, 0.99, 0.5, db, ax=ax,
                            asm1="asm1", asm2="asm2", pid_tie_eps=0.0001, binsize=1000)
    n = sum(len(c.get_segments()) for c in ax.collections if isinstance(c, LineCollection))
    plt.close(fig)
    check(n == 1, f"overlay should draw the clipped overhanging read, drew {n}")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_clip_tests_")
    tests = [
        (test_ref_to_read_inverts, False),
        (test_overlapping_segment_picks_largest, False),
        (test_visible_span_forward, False),
        (test_visible_span_reverse, False),
        (test_fully_contained_unchanged, False),
        (test_block_outside_still_dropped, False),
        (test_map_in_segment_clamps, False),
        (test_clipped_line_keeps_slope, True),
        (test_overlay_recovers_overhanging_read, True),
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
