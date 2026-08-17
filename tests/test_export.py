"""Tests for ViewState and figure export (Phase 4).

The risky part of export is the coordinate handover: qaview works in exact-bp plot
space, `report_plot` pads segments to whole bins.  `ViewState` bridges them by
carrying genomic intervals only, so these tests pin the clipping arithmetic and the
round trip rather than the rendering.

Run: venv/bin/python tests/test_export.py
"""

import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

import pysam

from fixtures import READS_A, READS_B, SEGS_A, SEGS_B, make_region_data, write_bam, write_db

from qatool.SegmentStorage import Segment
from qaview.state import ExportJob, ViewState, visible_segments

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


SEGS_INFO = [("ctgA", 1_000, 5_000, 0), ("ctgB", 200, 1_200, 4_000)]


def test_visible_segments_full():
    got = visible_segments(SEGS_INFO, "x")
    check(got == [Segment("x", "ctgA", 1_000, 5_000), Segment("x", "ctgB", 200, 1_200)],
          f"unbounded clip must reproduce the region exactly, got {got}")


def test_visible_segments_window():
    got = visible_segments(SEGS_INFO, "x", 500, 4_500)
    check(got == [Segment("x", "ctgA", 1_500, 5_000), Segment("x", "ctgB", 200, 700)],
          f"window must map back to genomic coords, got {got}")


def test_visible_segments_outside():
    check(visible_segments(SEGS_INFO, "x", 9_000, 9_999) == [], "window past the end -> nothing")
    check(visible_segments(SEGS_INFO, "x", -500, 0) == [], "window before the start -> nothing")


def test_visible_segments_rounds_not_truncates():
    """Bokeh hands back floats; truncating would shave a base off every export."""
    got = visible_segments(SEGS_INFO, "x", 0.6, 3_999.4)
    check(got[0].start == 1_001 and got[0].end == 4_999,
          f"expected rounded 1001..4999, got {got[0].start}..{got[0].end}")


def test_visible_segments_clamps_to_segment():
    got = visible_segments([("c", 100, 200, 0)], "x", -50, 5_000)
    check(got == [Segment("x", "c", 100, 200)],
          f"an over-wide window must not invent coordinates, got {got}")


def _state(tmp, rd, **kw):
    return ViewState.from_region(
        rd, db="db", bam1="b1", bam2="b2", fa1="f1", fa2="f2",
        read_names={"r1", "r2"}, pid_threshold=0.99, plen_threshold=0.9,
        best_pid_delta=0.001, pid_tie_eps=0.0001, **kw)


def test_from_region_full_and_clipped(tmp):
    rd = make_region_data(tmp, "vs")
    full = _state(tmp, rd)
    check(len(full.segments) == len(rd.segments),
          f"full view must keep every segment, got {len(full.segments)}")
    check({s.assembly for s in full.segments} == {rd.asm1, rd.asm2},
          "both assemblies must be present or report_plot cannot draw")

    half = _state(tmp, rd, x_range=(0, rd.plot_w / 2), y_range=(0, rd.plot_h / 2))
    span_full = sum(s.end - s.start for s in full.segments)
    span_half = sum(s.end - s.start for s in half.segments)
    check(span_half < span_full, f"clipped view must be smaller: {span_half} vs {span_full}")


def test_state_json_roundtrip(tmp):
    rd = make_region_data(tmp, "js")
    st = _state(tmp, rd, x_range=(10, rd.plot_w - 10))
    back = ViewState.from_json(st.to_json())
    check(back.segments == st.segments, "segments must survive the JSON round trip")
    check(all(isinstance(s, Segment) for s in back.segments),
          "segments must come back as Segment, not lists")
    check((back.asm1, back.asm2, back.dpi, back.binsize) ==
          (st.asm1, st.asm2, st.dpi, st.binsize), "scalar fields must round-trip")
    json.loads(st.to_json())          # must be valid JSON for the worker to read


def test_binsize_defaults_to_report_plot(tmp):
    """qaview's auto binsize has no k floor; report_plot's does, so let it decide."""
    rd = make_region_data(tmp, "bs", binsize=None)          # qaview's auto rule
    st = _state(tmp, rd)
    check(st.binsize is None,
          "export must default binsize to None so report_plot applies its max(k, ...) rule")
    # The hazard is concrete: qaview auto-sizes to ~1000 bins with no k floor, so on a
    # small region it lands below k and the k-mer mosaic would refuse to build.
    check(rd.binsize < st.k,
          f"expected auto binsize below k to demonstrate the hazard; got {rd.binsize} vs k={st.k}")


def test_describe(tmp):
    rd = make_region_data(tmp, "de")
    text = _state(tmp, rd).describe()
    check(rd.asm1 in text and rd.asm2 in text, f"describe must name both assemblies: {text}")
    check("bp" in text, f"describe must report spans: {text}")


def test_export_job_writes_state(tmp):
    rd = make_region_data(tmp, "ej")
    st = _state(tmp, rd)
    work = os.path.join(tmp, "exportwork")
    job = ExportJob(st, os.path.join(tmp, "out.png"), work, python="/nonexistent/python")
    check(job.argv[:2] == ["-m", "qaview.export_worker"], f"unexpected argv {job.argv}")
    check(job.label == "export", job.label)
    # start() writes the state file before launching; the bogus interpreter then fails.
    try:
        job.start()
    except FileNotFoundError:
        pass
    state_file = os.path.join(work, "export_state.json")
    check(os.path.exists(state_file), "state must be written for the child to read")
    if os.path.exists(state_file):
        with open(state_file) as fh:
            check(ViewState.from_json(fh.read()).segments == st.segments,
                  "written state must match")


def test_export_worker_renders(tmp):
    """Full path: state file in, figure out, through the real report_plot."""
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

    st = ViewState(db=db, bam1=bam1, bam2=bam2, fa1=fa, fa2=fa,
                   asm1="asm1", asm2="asm2",
                   segments=[SEGS_A[0], SEGS_B[0]], read_names=["r1"],
                   pid_threshold=0.99, plen_threshold=0.5,
                   best_pid_delta=0.001, pid_tie_eps=0.0001, dpi=60)
    sf = os.path.join(tmp, "st.json")
    with open(sf, "w") as fh:
        fh.write(st.to_json())
    out = os.path.join(tmp, "fig.png")

    r = subprocess.run([sys.executable, "-m", "qaview.export_worker", sf, out],
                       capture_output=True, text=True, timeout=600)
    check(r.returncode == 0, f"worker failed:\n{r.stdout[-800:]}\n{r.stderr[-800:]}")
    check(os.path.exists(out) and os.path.getsize(out) > 5_000,
          "export must produce a real figure")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_export_tests_")
    tests = [
        (test_visible_segments_full, False),
        (test_visible_segments_window, False),
        (test_visible_segments_outside, False),
        (test_visible_segments_rounds_not_truncates, False),
        (test_visible_segments_clamps_to_segment, False),
        (test_from_region_full_and_clipped, True),
        (test_state_json_roundtrip, True),
        (test_binsize_defaults_to_report_plot, True),
        (test_describe, True),
        (test_export_job_writes_state, True),
        (test_export_worker_renders, True),
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
