"""Structural tests for the qaview Bokeh layer (Phase 1).

Synthetic fixtures only.  These assert the properties the interactive behaviour rests
on — shared range objects, a shared crosshair instance, square dotplot pixels — since
those are what silently break and cannot be seen in a screenshot.

Run: venv/bin/python tests/test_glyphs.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

from bokeh.models import CrosshairTool, GridPlot, WheelZoomTool

from fixtures import make_region_data

from qaview.glyphs import (
    BAR_SIZE, MAIN_WIDTH, MARGIN_SIZE, contig_bar_figure, contig_colors,
    coverage_figure, dotplot_figure, main_height, report_layout, save_html,
)

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


def test_wheel_zoom_on_dotplot_only(tmp):
    rd = make_region_data(tmp, "wz")
    p = _panels(report_layout(rd))
    check(isinstance(p["main"].toolbar.active_scroll, WheelZoomTool),
          "wheel zoom must be active on the dotplot")
    for name in ("cov_top", "cov_right"):
        wheels = [t for t in p[name].tools if isinstance(t, WheelZoomTool)]
        check(not wheels, f"{name} must not have wheel zoom; it would desync the axes")


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
    # + 1 for the invisible per-bin hover target appended last.
    check(len(p.renderers) == n_series + 1,
          f"expected {n_series} stacked areas + 1 hit layer, got {len(p.renderers)}")
    check(abs(p.x_range.end - float(rd.cov1.total().max())) < 1e-9,
          f"value axis {p.x_range.end} should reach max total {rd.cov1.total().max()}")

    # Each bin contributes both its edges, so the area steps instead of interpolating.
    left, right = rd.bin_edges(rd.asm1)
    expected = [v for pair in zip(left, right) for v in pair]
    data = p.renderers[-2].data_source.data
    check(list(data["y"]) == expected,
          f"area must step on exact-bp bin edges; got {list(data['y'])[:6]}...")

    # Top of the last stacked area is the total coverage, doubled per bin.
    total = rd.cov1.total()
    check([round(v, 9) for v in data["x2"][0::2]] == [round(v, 9) for v in total],
          "top of the stack must equal summed coverage")


def test_coverage_has_no_bin_borders(tmp):
    """Bins must not be separated by drawn edges."""
    rd = make_region_data(tmp, "nb")
    for orient, asm in (("x", rd.asm2), ("y", rd.asm1)):
        p = coverage_figure(rd, asm, orient)
        for r in p.renderers[:-1]:                     # skip the invisible hit layer
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


def test_color_modes(tmp):
    rd = make_region_data(tmp, "cm")
    from qaview.glyphs import SUPPLEMENTARY_COLOR
    from qatool.plotters.colors import SUPPORT_COLORS

    from qaview.glyphs import LAYER_STYLE
    allowed = ({c for c, _ in LAYER_STYLE.values() if c is not None}
               | set(SUPPORT_COLORS.values()))
    by_class = dotplot_figure(rd, color_by="class").renderers[0].data_source.data["color"]
    check(set(by_class) <= allowed,
          f"class colours must come from the layer palette, got {set(by_class)}")

    by_block = dotplot_figure(rd, color_by="block").renderers[0].data_source.data["color"]
    check(set(by_block) <= {"black", SUPPLEMENTARY_COLOR},
          f"block mode is black + supplementary red, got {set(by_block)}")

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


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_glyph_tests_")
    tests = [test_square_pixels, test_ranges_are_linked, test_shared_crosshair,
             test_wheel_zoom_on_dotplot_only, test_panel_proportions,
             test_coverage_areas_stack, test_coverage_has_no_bin_borders,
             test_contig_colors_and_bars,
             test_color_modes, test_save_html]
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
