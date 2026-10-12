"""Structural tests for the qaview Bokeh layer (Phase 1).

Synthetic fixtures only.  These assert the properties the interactive behaviour rests
on — shared range objects, a shared crosshair instance, square dotplot pixels — since
those are what silently break and cannot be seen in a screenshot.

Run: venv/bin/python tests/test_glyphs.py
"""

import os
import sys

import numpy as np
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

from bokeh.models import CrosshairTool, GridPlot, WheelZoomTool

from fixtures import make_region_data

from qaview.glyphs import (
    BAR_SIZE, MAIN_WIDTH, MARGIN_SIZE, contig_bar_figure, contig_colors,
    END_MARK, contig_end_marks, coverage_figure, dotplot_figure, fit_main, main_height, report_layout, save_html,
    hit_renderer,
)
from qaview.region import marginal_bins

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _grid(layout):
    """The GridPlot inside the returned column (which also holds the readout Div)."""
    if isinstance(layout, GridPlot):
        return layout
    for child in layout.children:
        if isinstance(child, GridPlot):
            return child
    raise AssertionError("no GridPlot in layout")


def _panels(layout):
    """name -> figure, recovered from the gridplot by position."""
    by_pos = {(r, c): child for child, r, c in _grid(layout).children}
    return {"cov_top": by_pos[(0, 1)], "ybar": by_pos[(1, 0)],
            "main": by_pos[(1, 1)], "cov_right": by_pos[(1, 2)],
            "xbar": by_pos[(2, 1)]}


def test_square_pixels(tmp):
    """Dotplot aspect must track plot_h/plot_w or diagonals sit off 45 degrees."""
    rd = make_region_data(tmp, "sq")
    h = main_height(rd, MAIN_WIDTH)
    expected = round(MAIN_WIDTH * rd.plot_h / rd.plot_w)
    check(h == expected, f"height {h} != {expected} for {rd.plot_w}x{rd.plot_h}")
    check(h != MAIN_WIDTH, "fixture is deliberately non-square; a square figure would distort")

    p = dotplot_figure(rd)
    check(p.height == h, f"dotplot height {p.height} != {h}")
    check(p.x_range.end == rd.plot_w and p.y_range.end == rd.plot_h,
          "dotplot ranges must span the full region")


def test_fit_main_stays_in_its_box(tmp):
    """The dotplot never exceeds the box it is given, however small or lopsided.

    main_height clamps (200..4000 px), so recomputing the width from the height
    once overflowed a narrow box on a very wide locus; a 200 px floor overflowed
    small windows.
    """
    from types import SimpleNamespace

    for pw, ph in ((1_000_000, 1_000), (1_000, 1_000_000), (300_000, 200_000)):
        rd = SimpleNamespace(plot_w=pw, plot_h=ph)
        for bw, bh in ((1100, 800), (300, 150), (90, 70), (60, 400), (500, 60)):
            w, h = fit_main(rd, bw, bh)
            check(w <= bw and h <= bh, f"{pw}x{ph} locus in {bw}x{bh}: got {w}x{h}")


def test_ranges_are_linked(tmp):
    """Panels must share the *same* Range1d objects, not equal copies."""
    rd = make_region_data(tmp, "lk")
    p = _panels(report_layout(rd))
    check(p["cov_top"].x_range is p["main"].x_range, "top coverage x_range not linked")
    check(p["xbar"].x_range is p["main"].x_range, "bottom contig bar x_range not linked")
    check(p["cov_right"].y_range is p["main"].y_range, "right coverage y_range not linked")
    check(p["ybar"].y_range is p["main"].y_range, "left contig bar y_range not linked")
    # The marginals' value axes must stay independent of the dotplot's.
    check(p["cov_top"].y_range is not p["main"].y_range, "coverage y_range must be its own")
    check(p["cov_right"].x_range is not p["main"].x_range, "coverage x_range must be its own")


def test_shared_crosshair(tmp):
    """One CrosshairTool instance across panels is what makes the ruler project."""
    rd = make_region_data(tmp, "ch")
    p = _panels(report_layout(rd))
    tools = []
    for name in ("main", "cov_top", "cov_right"):
        found = [t for t in p[name].tools if isinstance(t, CrosshairTool)]
        check(len(found) == 1, f"{name} should have exactly one crosshair, got {len(found)}")
        tools.extend(found)
    check(len({id(t) for t in tools}) == 1,
          "all panels must share ONE CrosshairTool instance for a linked ruler")
    spans = tools[0].overlay
    check(isinstance(spans, tuple) and len(spans) == 2, f"expected 2 span overlays, got {spans}")
    check({s.dimension for s in spans} == {"width", "height"}, "need both ruler directions")


def test_wheel_zoom_on_marginals_zooms_one_axis(tmp):
    """The dotplot's wheel zooms both axes; a marginal's only its genome axis -- the
    one it shares with the dotplot -- never the shared coverage scale."""
    rd = make_region_data(tmp, "wz")
    p = _panels(report_layout(rd))
    w = p["main"].toolbar.active_scroll
    check(isinstance(w, WheelZoomTool) and w.dimensions == "both",
          "wheel zoom must be active on the dotplot, both axes")
    for name, dim, shared in (("cov_top", "width", "x_range"), ("cov_right", "height", "y_range")):
        wheels = [t for t in p[name].tools if isinstance(t, WheelZoomTool)]
        check(len(wheels) == 1 and wheels[0].dimensions == dim and p[name].toolbar.active_scroll is wheels[0],
              f"{name}: one active wheel zoom on its genome axis ({dim})")
        check(getattr(p[name], shared) is getattr(p["main"], shared),
              f"{name}: that axis is the dotplot's, so the dotplot follows")


def test_panel_proportions(tmp):
    """The flagged Phase 1 risk: gridplot must honour the report_plot proportions."""
    rd = make_region_data(tmp, "pp")
    layout = report_layout(rd)
    check(isinstance(_grid(layout), GridPlot), "layout must contain a GridPlot")
    p = _panels(layout)
    h = main_height(rd, MAIN_WIDTH)
    check(p["main"].width == MAIN_WIDTH and p["main"].height == h, "dotplot size")
    check(p["cov_top"].width == MAIN_WIDTH and p["cov_top"].height == MARGIN_SIZE,
          "top marginal must be as wide as the dotplot")
    check(p["cov_right"].height == h and p["cov_right"].width == MARGIN_SIZE,
          "right marginal must be as tall as the dotplot")
    check(p["xbar"].width == MAIN_WIDTH and p["xbar"].height == BAR_SIZE, "bottom contig bar")
    check(p["ybar"].height == h and p["ybar"].width == BAR_SIZE, "left contig bar")


def test_coverage_areas_stack(tmp):
    """Stacked areas must reach the summed coverage and step on exact-bp bin edges."""
    rd = make_region_data(tmp, "cq")
    p = coverage_figure(rd, rd.asm1, "y")
    n_series = len(rd.bucket_keys) + 1
    # + the dropout marks and the invisible per-bin hover target, appended last.
    check(len(p.renderers) == n_series + 2,
          f"expected {n_series} stacked areas + dropouts + hit layer, got {len(p.renderers)}")
    total = hit_renderer(p).data_source.data["total"]
    check(abs(p.x_range.end - float(total.max())) < 1e-9,
          f"value axis {p.x_range.end} should reach max total {total.max()}")

    # Each bin contributes both its edges, so the area steps instead of interpolating.
    b = marginal_bins(rd.asm1_segs_info, rd.placements[rd.asm1], [0] * len(rd.read_names), 1,
                      0, rd.plot_h)
    expected = [v for pair in zip(b["left"], b["right"]) for v in pair]
    data = p.renderers[-3].data_source.data
    check(list(data["e"]) == expected,
          f"area must step on exact-bp bin edges; got {list(data['e'])[:6]}...")

    # Top of the last stacked area is the total coverage, doubled per bin.
    check([round(v, 9) for v in data["hi"][0::2]] == [round(v, 9) for v in b["series"][0]],
          "top of the stack must equal summed coverage")


def test_coverage_has_no_bin_borders(tmp):
    """Bins must not be separated by drawn edges."""
    rd = make_region_data(tmp, "nb")
    for orient, asm in (("x", rd.asm2), ("y", rd.asm1)):
        p = coverage_figure(rd, asm, orient)
        for r in p.renderers[:-2]:                     # skip the dropouts and hit layer
            g = r.glyph
            check(getattr(g, "line_alpha", 0) in (0, None) or
                  getattr(g, "line_color", None) is None,
                  f"{orient}: {type(g).__name__} must draw no border between bins")


def test_contig_colors_and_bars(tmp):
    rd = make_region_data(tmp, "cc")
    colors = contig_colors(rd)
    check(set(colors) == {"ctgA", "ctgB"}, f"expected both contigs, got {set(colors)}")
    check(all(c.startswith("#") for c in colors.values()), f"need hex colours, got {colors}")
    check(colors["ctgA"] != colors["ctgB"], "contigs must be distinguishable")

    bar = contig_bar_figure(rd, rd.asm1, "y", colors)
    check(len(bar.renderers) == 1, "asm1's two segments are both ctgA -> one merged run")
    check(bar.y_range.end == rd.plot_h, "contig bar must span the full axis")


def test_contig_end_marks(tmp):
    """A triangle where a segment reaches its contig's start or end, pointing inward."""
    import math

    rd = make_region_data(tmp, "ce")
    check(rd.contig_lengths == {"asm1": {"ctgA": 20_000}, "asm2": {"ctgB": 12_000}},
          f"contig lengths come from the BAM headers: {rd.contig_lengths}")
    check(contig_end_marks(rd, "asm1", tol=0) == [] and contig_end_marks(rd, "asm2", tol=0) == [],
          "the fixture's segments reach no contig end")
    near = contig_end_marks(rd, "asm2", tol=1_000)          # ctgB 1,000-4,300 of 12,000
    check([(m[3], m[5]) for m in near] == [("start", 1_000)],
          f"an edge within the tolerance counts, with how far short it is: {near}")

    # asm1: ctgA 2,000-6,500 and 9,000-11,200 -> as if ctgA started at 2,000 and
    # ended at 11,200 (the length says where the end is; the start is 0 by definition)
    rd.contig_lengths = {"asm1": {"ctgA": 11_200}, "asm2": {"ctgB": 4_300}}
    rd.asm1_segs_info = [("ctgA", 0, 6_500, 0)] + [tuple(x) for x in rd.asm1_segs_info[1:]]
    marks = contig_end_marks(rd, "asm1", tol=0)
    off2 = rd.asm1_segs_info[1][3]
    check([(m[0], m[1], m[3]) for m in marks] == [(0, 1, "start"), (off2 + 2_200, -1, "end")],
          f"asm1: start of the first segment, end of the second: {marks}")
    check([(m[1], m[3]) for m in contig_end_marks(rd, "asm2", tol=0)] == [(-1, "end")],
          "asm2: ctgB's end only")

    colors = contig_colors(rd)
    bar = contig_bar_figure(rd, rd.asm1, "y", colors)
    tri = [r for r in bar.renderers if getattr(r.glyph, "marker", None) == "triangle"]
    check(len(tri) == 1, f"one triangle renderer: {bar.renderers}")
    d = tri[0].data_source.data
    check(d["angle"] == [0.0, math.pi], f"y bar: up at the start, down at the end: {d['angle']}")
    check(d["pos"][0] > 0 and d["pos"][1] < off2 + 2_200,
          f"each triangle sits inside its segment: {d['pos']}")
    per_px = rd.plot_h / bar.height
    check(abs((d["pos"][0] - 0) / per_px - END_MARK / 2) < 1e-6, "by half its size, in pixels")
    xbar = contig_bar_figure(rd, rd.asm2, "x", colors)
    xd = [r for r in xbar.renderers if getattr(r.glyph, "marker", None) == "triangle"][0].data_source.data
    # ctgB 1,000-4,300 of 4,300: its end, and its start within the default 1 kb
    check(xd["angle"] == [-math.pi / 2, math.pi / 2],
          f"x bar: a start points right, an end left: {xd['angle']}")
    check(bar.y_range.js_property_callbacks.get("change:start"), "the offset follows zoom")


def test_color_modes(tmp):
    rd = make_region_data(tmp, "cm")
    from qaview.glyphs import SUPPLEMENTARY_COLOR
    from qatool.plotters.colors import SUPPORT_COLORS

    from qaview.glyphs import LAYER_STYLE, pack_rgba
    # colours travel packed (0xRRGGBBAA), one binary array
    allowed = {pack_rgba(c) for c in ({c for c, _ in LAYER_STYLE.values() if c is not None}
                                      | set(SUPPORT_COLORS.values()))}
    by_class = dotplot_figure(rd, color_by="class").renderers[0].data_source.data["color"]
    check(by_class.dtype == np.uint32 and set(by_class.tolist()) <= allowed,
          f"class colours must come from the layer palette, got {set(by_class.tolist())}")

    by_block = dotplot_figure(rd, color_by="block").renderers[0].data_source.data["color"]
    check(set(by_block.tolist()) <= {pack_rgba("black"), pack_rgba(SUPPLEMENTARY_COLOR)},
          f"block mode is black + supplementary red, got {set(by_block.tolist())}")

    try:
        dotplot_figure(rd, color_by="nonsense")
        check(False, "an unknown color_by must raise")
    except ValueError:
        pass


def test_save_html(tmp):
    rd = make_region_data(tmp, "sv")
    out = os.path.join(tmp, "report.html")
    save_html(report_layout(rd), out)
    check(os.path.exists(out), "no HTML written")
    body = open(out, encoding="utf-8").read()
    check(len(body) > 100_000, f"bokeh JS should be inlined; file is only {len(body)} bytes")
    check("<script" in body and "cdn.bokeh.org" not in body,
          "resources must be inline so the file works offline")


def test_marginals_share_one_scale(tmp):
    """As in report_plot: both coverage marginals are the same size on the same
    scale, so equal coverage is equally long on either."""
    rd = make_region_data(tmp, "ms")
    p = _panels(report_layout(rd))
    top, right = p["cov_top"], p["cov_right"]
    peak = max(float(np.max(hit_renderer(f).data_source.data["total"])) for f in (top, right))
    check(top.y_range is right.x_range, "one range object: zoom one, the other follows")
    check(abs(top.y_range.end - peak) < 1e-9 and top.y_range.start == 0,
          f"from 0 to the larger peak: {top.y_range.end} vs {peak}")
    check(top.height == right.width, f"same size: {top.height} vs {right.width}")
    for name in ("cov_top", "cov_right"):
        src = p[name].renderers[-1].data_source.data
        check(np.allclose(src["top"], peak), f"{name}: the hover target reaches the shared top")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_glyph_tests_")
    tests = [test_square_pixels, test_fit_main_stays_in_its_box, test_ranges_are_linked, test_shared_crosshair,
             test_wheel_zoom_on_marginals_zooms_one_axis, test_panel_proportions,
             test_coverage_areas_stack, test_coverage_has_no_bin_borders,
             test_contig_colors_and_bars, test_contig_end_marks,
             test_color_modes, test_save_html, test_marginals_share_one_scale]
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
