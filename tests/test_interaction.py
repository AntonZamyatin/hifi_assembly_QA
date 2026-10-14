"""Structural tests for the qaview interaction layer (Phase 2).

Hover and highlight run entirely in the browser, so the parts that can break silently
are the two CSR indexes and the read-name remapping between them.  Those are checked
here in Python, against the same arrays the CustomJS callbacks receive.

Run: venv/bin/python tests/test_interaction.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")

from bokeh.models import CustomJS, Div, GridPlot, HoverTool, TapTool

from fixtures import make_region_data

from qaview.glyphs import (
    coverage_figure, dotplot_figure, hit_renderer, line_renderer, report_layout,
)

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _panels(layout):
    grid = next(c for c in layout.children if isinstance(c, GridPlot))
    by_pos = {(r, c): child for child, r, c in grid.children}
    return {"cov_top": by_pos[(0, 1)], "main": by_pos[(1, 1)], "cov_right": by_pos[(1, 2)]}


def test_hit_layer_present_and_invisible(tmp):
    rd = make_region_data(tmp, "hit")
    for orient, asm in (("x", rd.asm2), ("y", rd.asm1)):
        p = coverage_figure(rd, asm, orient)
        hit = hit_renderer(p)
        n_bins = len(p.renderers[0].data_source.data["e"]) // 2
        check(len(hit.data_source.data["left"]) == n_bins,
              f"{orient}: hit layer must have one quad per bin")
        check(hit.glyph.fill_alpha == 0.0, f"{orient}: hit layer must be invisible")
        check(hit.selection_glyph is not None and hit.selection_glyph.fill_alpha > 0,
              f"{orient}: selected bins must become visible")
        data = hit.data_source.data
        check(all(c and e > s0 for c, s0, e in zip(data["contig"], data["bin_start"], data["bin_end"])),
              f"{orient}: hit layer bins must carry genomic identity for the tooltip")


def test_hover_tools_wired(tmp):
    rd = make_region_data(tmp, "hv")
    p = _panels(report_layout(rd))

    main_hovers = [t for t in p["main"].tools if isinstance(t, HoverTool)]
    check(len(main_hovers) == 1, f"dotplot needs one HoverTool, got {len(main_hovers)}")
    fields = {f for _, f in main_hovers[0].tooltips}
    check(any("@read{name}" in f for f in fields), "read tooltip must name the read")
    check(set(main_hovers[0].formatters) >= {"@read", "@cls", "@c1", "@c2", "@supp"},
          "the codes the rows hold are turned back into names")
    check(any("pid1" in f for f in fields) and any("pid2" in f for f in fields),
          "read tooltip must show pid on both assemblies")
    check(main_hovers[0] not in (p["main"].toolbar.active_inspect or []),
          "the read tooltip starts switched off (it covered the shift+drag rectangle)")

    for name in ("cov_top", "cov_right"):
        hovers = [t for t in p[name].tools if isinstance(t, HoverTool)]
        check(len(hovers) == 1, f"{name} needs one HoverTool, got {len(hovers)}")
        check(isinstance(hovers[0].callback, CustomJS),
              f"{name} hover must carry the bin->read callback")
        args = hovers[0].callback.args
        for key in ("lines", "M", "O", "reads"):
            check(key in args, f"{name} callback missing arg {key}")


def test_selection_glyphs(tmp):
    rd = make_region_data(tmp, "sel")
    p = _panels(report_layout(rd))
    r = line_renderer(p["main"])
    check(r.selection_glyph is not None, "selected reads need a distinct glyph")
    check(r.nonselection_glyph is not None, "unselected reads must dim")
    # Alpha is a per-row field so a layer hidden by a Tier-3 toggle stays hidden even
    # when selected; only the nonselection state is a hard zero.
    check(r.glyph.line_alpha == "alpha" and r.selection_glyph.line_alpha == "sel_alpha",
          "base and selected glyphs must honour the per-row alpha fields")
    src = r.data_source.data
    check(all(float(a) == float(s) for a, s, lay in zip(src["alpha"], src["sel_alpha"], src["layer"])
              if lay != 4), "an ordinary row's selected alpha is its alpha")
    check(r.nonselection_glyph.line_alpha == 0.0,
          "unselected reads must be fully hidden")
    check(any(isinstance(t, TapTool) for t in p["main"].tools),
          "tapping a read is how a selection is made by hand")


def test_reverse_link_registered(tmp):
    rd = make_region_data(tmp, "rev")
    p = _panels(report_layout(rd))
    src = line_renderer(p["main"]).data_source
    cbs = src.selected.js_property_callbacks.get("change:indices", [])
    check(len(cbs) == 1, f"expected one read->bin callback, got {len(cbs)}")
    if cbs:
        args = cbs[0].args
        check(len(args["margs"]) == 2,
              "read->bin callback must reach both coverage panels")


def test_draw_order_and_widths(tmp):
    """Rows must be stored back-to-front, so one glyph renderer paints the layers."""
    from qaview.glyphs import LAYER_STYLE
    from qaview.region import (
        LAYER_COMMON, LAYER_OTHER, LAYER_SUPPLEMENTARY, LAYER_SUPPORT, _layer_of,
    )

    check(LAYER_COMMON < LAYER_OTHER < LAYER_SUPPORT < LAYER_SUPPLEMENTARY,
          "layer order must be common < other < support < supplementary")
    check(_layer_of("common", False) == LAYER_COMMON, "common reads sit at the bottom")
    check(_layer_of("other", False) == LAYER_OTHER, "unclassified above common")
    check(_layer_of("asm1_best", False) == LAYER_SUPPORT, "support above unclassified")
    check(_layer_of("common", True) == LAYER_SUPPLEMENTARY,
          "supplementary must win over classification")

    widths = [LAYER_STYLE[i][1] for i in (LAYER_COMMON, LAYER_OTHER, LAYER_SUPPORT)]
    check(widths[0] > widths[1] >= widths[2],
          f"common must be the widest line, got {widths}")
    check(LAYER_STYLE[LAYER_SUPPORT][0] is None,
          "support lines take their colour from the read class, not the layer")

    rd = make_region_data(tmp, "ord")
    layers = list(rd.lines.layer)
    check(layers == sorted(layers), f"stored rows must be in draw order, got {layers}")

    src = dotplot_figure(rd).renderers[0].data_source.data
    check(len(src["width"]) == len(rd.lines), "every line needs its own width")
    for lay, w in zip(layers, src["width"]):
        check(w == LAYER_STYLE[int(lay)][1], f"layer {lay} should draw at {LAYER_STYLE[int(lay)][1]}")


def test_hover_links_both_panels(tmp):
    """Hovering one coverage panel must also drive the opposite one."""
    rd = make_region_data(tmp, "both")
    p = _panels(report_layout(rd))
    for name, partner in (("cov_top", "cov_right"), ("cov_right", "cov_top")):
        hover = next(t for t in p[name].tools if isinstance(t, HoverTool))
        args = hover.callback.args
        check(args["O"]["hit"] is hit_renderer(p[partner]).data_source,
              f"{name} must target {partner}'s hit layer")
        check(args["M"]["hit"] is hit_renderer(p[name]).data_source,
              f"{name} must stripe its own hovered bin")


def test_nonselected_reads_hidden(tmp):
    rd = make_region_data(tmp, "hid")
    r = line_renderer(_panels(report_layout(rd))["main"])
    check(r.nonselection_glyph.line_alpha == 0.0,
          "reads outside the hovered bin must be invisible, not merely dimmed")


def test_band_detection_data(tmp):
    """Hovering a stacked band narrows the highlight to that read class.

    The band is resolved client-side from the cursor's value-axis position against
    this bin's cumulative band tops, so those tops must be on the hit source and the
    read -> band map must agree with the coverage stacking order.
    """
    from qaview.glyphs import _read_bucket_index

    rd = make_region_data(tmp, "band")
    n_series = len(rd.bucket_keys) + 1

    for orient, asm in (("x", rd.asm2), ("y", rd.asm1)):
        p = coverage_figure(rd, asm, orient)
        data = hit_renderer(p).data_source.data
        n_bins = len(data["left"])
        for s in range(n_series):
            check(f"cum_class_{s}" in data, f"{orient}: missing band top cum_class_{s}")
            check(len(data[f"cum_class_{s}"]) == n_bins,
                  f"{orient}: cum_class_{s} must have one value per bin")
        # Bands must be non-decreasing and end at the total.
        for i in range(n_bins):
            tops = [data[f"cum_class_{s}"][i] for s in range(n_series)]
            check(tops == sorted(tops), f"{orient}: bin {i} band tops not monotonic: {tops}")
        check([round(v, 9) for v in data[f"cum_class_{n_series - 1}"]] ==
              [round(v, 9) for v in data["total"]],
              f"{orient}: last band must close at the total coverage")

    band = _read_bucket_index(rd)
    check(len(band) == len(rd.read_names), "one band per global read")
    check(all(0 <= b <= len(rd.bucket_keys) for b in band), f"band index out of range: {set(band)}")
    for i, name in enumerate(rd.read_names):
        cls = rd.read_class.get(name)
        if cls in rd.bucket_keys:
            check(band[i] == rd.bucket_keys.index(cls),
                  f"{name}: band {band[i]} != stack position of class {cls}")

    # The stacking order in coverage_figure is bucket_keys then the leftover series,
    # which is exactly what the band index assumes.
    p = _panels(report_layout(rd))
    args = next(t for t in p["cov_top"].tools if isinstance(t, HoverTool)).callback.args
    check(args["names"]["names_class"] == list(rd.bucket_keys) + ["unassigned"],
          f"band names must mirror the stack order, got {args['names']}")
    check(list(args["reads"].data["band_class"]) == band, "the read -> band map the hover uses")
    check(args["mode"].text == "class", "unscored: the support classes are on show")
    check(args["vertical"] is True, "top marginal stacks along y, so the cursor's y picks the band")
    args_r = next(t for t in p["cov_right"].tools if isinstance(t, HoverTool)).callback.args
    check(args_r["vertical"] is False, "right marginal stacks along x")


def test_common_colour_matches_marginals(tmp):
    from qaview.glyphs import LAYER_STYLE
    from qaview.region import LAYER_COMMON
    from qatool.plotters.colors import SUPPORT_COLORS
    check(LAYER_STYLE[LAYER_COMMON][0] == SUPPORT_COLORS["common"],
          "common reads must use the same colour as the marginal stack")


def test_coordinate_readout(tmp):
    rd = make_region_data(tmp, "co")
    layout = report_layout(rd)
    divs = [c for c in layout.children if isinstance(c, Div) and c.visible]
    check(len(divs) == 2, f"expected coordinate + bin-status Divs, got {len(divs)}")
    p = _panels(layout)
    events = p["main"].js_event_callbacks.get("mousemove", [])
    check(len(events) == 1, "dotplot must emit a mousemove for the coordinate readout")
    if events:
        args = events[0].args
        check("ysegs" in args and "xsegs" in args,
              "readout needs both axis maps to convert plot -> genome")
        # The JS lookup must mirror plot_to_genome, including its segment bounds.
        check([list(s) for s in args["ysegs"]] ==
              [[c, s, e, o, False] for c, s, e, o in rd.asm1_segs_info],
              "ysegs must match asm1_segs_info exactly")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_interaction_tests_")
    tests = [test_hit_layer_present_and_invisible,
             test_hover_tools_wired, test_selection_glyphs, test_reverse_link_registered,
             test_draw_order_and_widths, test_hover_links_both_panels,
             test_nonselected_reads_hidden,
             test_band_detection_data, test_common_colour_matches_marginals,
             test_coordinate_readout]
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
