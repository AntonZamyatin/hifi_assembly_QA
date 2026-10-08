"""Tests for the pair-HMM score in qaview: log-ratio colouring and the Delta S panel.

The colouring has one hard constraint, inherited from the hover and selection
callbacks: rows are addressed by position, so a colour mode may change colours,
widths and opacities but never the number or order of rows. Most checks here
are about that, plus the promise that the clamp C is display-only.

Run: venv/bin/python tests/test_score_view.py
"""

import math
import os
import pickle
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

from fixtures import make_region_data, write_bam

import numpy as np

from qatool.plotters.coverage import compute_coverage_bins
from qatool.SegmentStorage import Segment

from qamodel.exposure import Exposure
from qamodel.score import LocusScore, ReadScore
from qaview.app import QaviewApp, dataset_paths
from qaview.glyphs import (UNSCORED_COLOR, _line_source, class_stacking, ratio_band,
                           ratio_color, ratio_stacking, ratio_strength, report_layout)
from qaview.score import ScoreJob, score_cache_path

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def fake_score(names, ratios, dropped=()):
    reads = [ReadScore(n, 1000, 0.0, -r, 1, 1, 0.0) for n, r in zip(names, ratios)]
    x = Exposure(1000, 1000, 1, 1, 500.0, 500.0, 1e-3)
    return LocusScore("asm1", "asm2", reads, x, C=200.0, dropped=list(dropped))


def test_strength_scale():
    check(ratio_strength(0.0, 200) == 0.0, "zero must be neutral")
    check(abs(ratio_strength(5, 200) + ratio_strength(-5, 200)) < 1e-12, "must be symmetric")
    check(ratio_strength(math.inf, 200) == 1.0 and ratio_strength(-math.inf, 200) == -1.0,
          "one-sided reads must saturate")
    check(ratio_strength(10_000, 200) == 1.0, "must saturate at C")
    check(ratio_strength(float("nan"), 200) == 0.0, "NaN must not colour")
    ts = [ratio_strength(r, 200) for r in (0.1, 1, 10, 100)]
    check(ts == sorted(ts), "must be monotone")
    # the point of the log scale: a decisive -27 must not look neutral at C = 200
    check(abs(ratio_strength(-27, 200)) > 0.5,
          f"-27 nats should read as clearly coloured, got {ratio_strength(-27, 200):.2f}")
    check(ratio_color(0.0) == "#d2d2ca", "zero maps to the neutral grey")


def _rgb(packed):
    return (int(packed) >> 24) & 255, (int(packed) >> 16) & 255, (int(packed) >> 8) & 255


def test_line_source_keeps_rows(tmp):
    from qaview.glyphs import _line_sources, pack_rgba

    rd = make_region_data(tmp, "lr")
    base = _line_source(rd, "class")
    names = sorted(set(rd.read_names))
    ratios = {names[0]: -27.0, names[1]: 3.0}
    src, reads, _ = _line_sources(rd, "log-ratio", ratios=ratios, dropped={names[2]}, C=200)
    n = len(base.data["x0"])
    check(len(src.data["x0"]) == n and list(src.data["read"]) == list(base.data["read"]),
          "log-ratio mode must keep every row in the same order")
    for g, colour in zip(src.data["read"], src.data["color"]):
        name, txt = reads.data["name"][g], reads.data["ratio_txt"][g]
        if name == names[0]:
            r, gr, b = _rgb(colour)
            check(gr > r, f"a read favouring asm2 must lean green, got {_rgb(colour)}")
            check(txt == "-27.00 nats", f"hover text {txt!r}")
        elif name == names[1]:
            r, gr, b = _rgb(colour)
            check(r > gr, f"a read favouring asm1 must lean red, got {_rgb(colour)}")
        elif name == names[2]:
            check(colour == pack_rgba(UNSCORED_COLOR) and txt.startswith("dropped"),
                  f"a dropped read must be grey and say so: {_rgb(colour)} {txt!r}")


def test_packed_ratio_colours_match_ratio_color():
    """The array path (and the page, which mirrors it) colours as ratio_color does."""
    import numpy as np

    from qaview.glyphs import pack_rgba, ratio_rgba, ratio_strengths

    rs = [-math.inf, -1e4, -27.0, -3.0, -0.4, 0.0, 0.4, 3.0, 27.0, 1e4, math.inf, float("nan")]
    for C in (50, 200, 1000):
        ts = ratio_strengths(np.array(rs), C)
        check(all(abs(t - ratio_strength(r, C)) < 1e-12 for t, r in zip(ts, rs)),
              f"strengths at C={C}: {ts}")
        got = ratio_rgba(ts)
        check(all(int(p) == pack_rgba(ratio_color(ratio_strength(r, C))) for p, r in zip(got, rs)),
              f"colours at C={C}")


def test_hidden_layers_stay_hidden(tmp):
    rd = make_region_data(tmp, "hl")
    ratios = {n: -10.0 for n in rd.read_names}
    src = _line_source(rd, "log-ratio", visible_layers=set(), ratios=ratios)
    check(all(float(a) == 0.0 for a in src.data["alpha"]),
          "hiding every layer must hide every line in log-ratio mode too")


def test_unscored_before_scoring(tmp):
    rd = make_region_data(tmp, "us")
    from qaview.glyphs import _line_sources, pack_rgba

    src, reads, _ = _line_sources(rd, "log-ratio", ratios=None)
    check(all(int(c) == pack_rgba(UNSCORED_COLOR) for c in src.data["color"]),
          "before scoring finishes every line is the unscored grey")
    check(all(t == "not scored yet" for t in reads.data["ratio_txt"]), "hover says so")
    report_layout(rd, color_by="log-ratio")                 # must not raise


def test_score_cache_ignores_display(tmp):
    p = dict(dataset_paths("b", "l", "asm1", "asm2"), asm1="asm1", asm2="asm2",
             contig="c", start=0, end=100, pid_threshold=0.997, plen_threshold=0.999,
             best_pid_delta=0.001, merge_gap=1, extension_window=1, ext_max=1,
             ext_rounds=1, max_rounds=1)
    a = score_cache_path(tmp, p)
    check(a == score_cache_path(tmp, dict(p, binsize=5000, color="log-ratio")),
          "display settings must not invalidate a score")
    check(a != score_cache_path(tmp, dict(p, asm2="other")),
          "a different assembly pair must be a different score")


def test_score_job_reads_cache(tmp):
    p = dict(dataset_paths("b", "l", "asm1", "asm2"), asm1="asm1", asm2="asm2",
             contig="c", start=0, end=100, pid_threshold=0.997, plen_threshold=0.999,
             best_pid_delta=0.001, merge_gap=1, extension_window=1, ext_max=1,
             ext_rounds=1, max_rounds=1)
    cache = os.path.join(tmp, "sc")
    os.makedirs(cache, exist_ok=True)
    ls = fake_score(["x"], [-5.0])
    with open(score_cache_path(cache, p), "wb") as fh:
        pickle.dump(ls, fh)
    job = ScoreJob(p, "unused.pkl", cache).start()
    check(job.done and job.from_cache and job.result.reads[0].ratio == -5.0,
          "a cached score must load without starting a subprocess")


def test_clamp_is_display_only(tmp):
    """Moving C changes Delta S in the panel and the colours, never the stored ratios."""
    app = QaviewApp(cache_dir=os.path.join(tmp, "cache"))
    app.rd = make_region_data(tmp, "cl")
    names = sorted(set(app.rd.read_names))
    app.score = fake_score(names, [-150.0] + [0.1] * (len(names) - 1))
    app.color_by.value = "log-ratio"
    app.clamp.value = 200
    app._render()
    html200 = app.score_pane.object
    app.clamp.value = 50
    app._render()
    html50 = app.score_pane.object
    want200 = f"{app.score.at(200):+.2f}"
    want50 = f"{app.score.at(50):+.2f}"
    check(want200 in html200 and want50 in html50 and want50 != want200,
          f"panel must show Delta S at the slider's C: {want200} / {want50}")
    check(app.score.reads[0].ratio == -150.0, "the slider must not touch stored ratios")


def test_rebucket_counts_repeated_placements(tmp):
    """Two tied alignments of one read in the same bin contribute 2/n there, not 1/n.

    This is why the bins carry per-read weights: regrouping from membership alone
    would halve such a read's coverage in exactly the repeats this tool is for.
    """
    bam = write_bam(os.path.join(tmp, "rep.bam"), [
        ("m", "ctgA", 2_100, 200, 0, 0),
        ("m", "ctgA", 2_400, 200, 0, 256),          # same bin, equally good
        ("u", "ctgA", 3_100, 300, 0, 0),
    ])
    cov = compute_coverage_bins(bam, [Segment("asm1", "ctgA", 2_000, 4_000)],
                                pid_threshold=0.99, plen_threshold=0.5, binsize=1000)
    m = cov.read_names.index("m")
    entry = [k for k in range(cov.bin_offsets[0], cov.bin_offsets[1])
             if cov.read_idx[k] == m]
    check(len(entry) == 1 and abs(cov.read_w[entry[0]] - 1.0) < 1e-12,
          f"read m should weigh 2 x 1/2 = 1.0 in bin 0, got {cov.read_w[entry]}")
    bands = cov.rebucket([0] * len(cov.read_names), 1)
    check(np.allclose(bands.sum(axis=0), cov.total()),
          f"rebucketed total {bands.sum(axis=0)} != {cov.total()}")


def test_ratio_bands():
    check([ratio_band(x) for x in (None, float("nan"), 0.5, -0.99, -1, -9.9, -10, -1e9,
                                   1, 12, math.inf)]
          == [5, 5, 0, 0, 1, 1, 2, 2, 3, 4, 4],
          "bands: indifferent < 1, weak 1-10, strong >= 10, per side")


def test_ratio_stacking_is_exact(tmp):
    rd = make_region_data(tmp, "rs")
    names = sorted(set(rd.read_names))
    ratios = {names[0]: -27.0, names[1]: 3.0}           # every other read unscored
    rs, cs = ratio_stacking(rd, ratios), class_stacking(rd)
    check(len(rs.keys) == len(rs.colors) == 6, "six bands, six colours")
    for asm in (rd.asm1, rd.asm2):
        a = np.sum(rs.series[asm], axis=0)
        b = np.sum(np.array(cs.series[asm], dtype=float), axis=0)
        check(np.allclose(a, b),
              f"{asm}: log-ratio bands must sum to the same coverage as the classes")
    got = [rs.read_band[rd.read_names.index(n)] for n in names]
    check(got == [2, 3] + [5] * (len(names) - 2),
          f"-27 is strong for asm2, +3 weak for asm1, the rest unscored; got {got}")


def test_layout_marginals_follow_colour(tmp):
    rd = make_region_data(tmp, "lm")
    names = sorted(set(rd.read_names))
    ratios = {n: -12.0 for n in names}

    def shown(layout):
        """The band names of the stacking on show, as the hover reads them."""
        from bokeh.models import CustomJS, HoverTool
        out = set()
        for m in layout.select({"type": HoverTool}):
            if isinstance(m.callback, CustomJS) and "names" in m.callback.args:
                a = m.callback.args
                out.add(tuple(a["names"][f"names_{a['mode'].text}"]))
        return out

    def drawn(layout):
        """Stackings drawn, and which are visible (the page switches them)."""
        from bokeh.models import GlyphRenderer
        tags = {}
        for r in layout.select({"type": GlyphRenderer}):
            if r.tags and str(r.tags[0]).startswith("stack:"):
                tags.setdefault(r.tags[0], set()).add(r.visible)
        return tags

    lay = report_layout(rd, color_by="log-ratio", ratios=ratios)
    lr = shown(lay)
    check(len(lr) == 1 and "not scored" in next(iter(lr)),
          f"log-ratio colouring must stack the marginals by log-ratio band, got {lr}")
    check(drawn(lay) == {"stack:ratio": {True}, "stack:class": {False}},
          f"both stackings drawn, the ratio one shown: {drawn(lay)}")
    lay = report_layout(rd, color_by="class", ratios=ratios)
    cl = shown(lay)
    check(len(cl) == 1 and "unassigned" in next(iter(cl)),
          f"class colouring must keep the class stacking, got {cl}")
    check(drawn(lay) == {"stack:ratio": {False}, "stack:class": {True}}, f"{drawn(lay)}")
    lay = report_layout(rd, color_by="log-ratio", ratios=None)
    pre = shown(lay)
    check(len(pre) == 1 and "unassigned" in next(iter(pre)),
          "before scoring, log-ratio mode falls back to class marginals")
    check(drawn(lay) == {"stack:class": {True}}, f"no ratio stacking before scoring: {drawn(lay)}")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_score_view_")
    tests = [(test_strength_scale, False), (test_line_source_keeps_rows, True),
             (test_hidden_layers_stay_hidden, True), (test_unscored_before_scoring, True),
             (test_score_cache_ignores_display, True), (test_score_job_reads_cache, True),
             (test_clamp_is_display_only, True),
             (test_rebucket_counts_repeated_placements, True), (test_ratio_bands, False),
             (test_ratio_stacking_is_exact, True), (test_layout_marginals_follow_colour, True),
             (test_packed_ratio_colours_match_ratio_color, False)]
    for fn, needs in tests:
        print(f"{fn.__name__} ...")
        fn(tmp) if needs else fn()
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
