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


def test_read_to_line_csr(tmp):
    """ReadLines.by_read must invert read_idx exactly."""
    rd = make_region_data(tmp, "csr")
    off, idx = rd.lines.by_read(len(rd.read_names))
    check(len(off) == len(rd.read_names) + 1, "offsets must be n_reads + 1")
    check(int(off[-1]) == len(rd.lines), "offsets must end at the line count")
    for r in range(len(rd.read_names)):
        rows = set(idx[off[r]:off[r + 1]].tolist())
        expected = {i for i, ri in enumerate(rd.lines.read_idx) if ri == r}
        check(rows == expected,
              f"read {rd.read_names[r]}: CSR rows {rows} != {expected}")


def test_bin_to_line_roundtrip(tmp):
    """Reproduce the CustomJS lookup in Python: bin -> reads -> line rows."""
    rd = make_region_data(tmp, "rt")
    cov = rd.cov1
    off, idx = rd.lines.by_read(len(rd.read_names))
    pos = {n: i for i, n in enumerate(rd.read_names)}

    checked = 0
    for b in range(len(cov.bins)):
        names = cov.reads_in_bin(b)
        if not names:
            continue
        rows = set()
        for n in names:
            g = pos.get(n)
            if g is None:
                continue
            rows.update(idx[off[g]:off[g + 1]].tolist())
        # Every selected line must belong to a read that really is in this bin.
        for r in rows:
            check(rd.read_names[rd.lines.read_idx[r]] in names,
                  f"bin {b} selected a line whose read is not in the bin")
        checked += 1
    check(checked > 0, "fixture produced no populated bins to test")


def test_coverage_read_names_are_a_subset(tmp):
    """The JS remap assumes every coverage read name exists in rd.read_names."""
    rd = make_region_data(tmp, "sub")
    for label, cov in (("cov1", rd.cov1), ("cov2", rd.cov2)):
        missing = set(cov.read_names) - set(rd.read_names)
        check(not missing,
              f"{label} has reads absent from RegionData.read_names: {missing}")


def test_hit_layer_present_and_invisible(tmp):
    rd = make_region_data(tmp, "hit")
    for orient, asm in (("x", rd.asm2), ("y", rd.asm1)):
        p = coverage_figure(rd, asm, orient)
        hit = hit_renderer(p)
        cov = rd.cov1 if asm == rd.asm1 else rd.cov2
        check(len(hit.data_source.data["bin_index"]) == len(cov.bins),
              f"{orient}: hit layer must have one quad per bin")
        check(hit.glyph.fill_alpha == 0.0, f"{orient}: hit layer must be invisible")
        check(hit.selection_glyph is not None and hit.selection_glyph.fill_alpha > 0,
              f"{orient}: selected bins must become visible")
        check(list(hit.data_source.data["contig"]) == list(cov.bin_contig),
              f"{orient}: hit layer bins must carry genomic identity for the tooltip")


def test_hover_tools_wired(tmp):
    rd = make_region_data(tmp, "hv")
    p = _panels(report_layout(rd))

    main_hovers = [t for t in p["main"].tools if isinstance(t, HoverTool)]
    check(len(main_hovers) == 1, f"dotplot needs one HoverTool, got {len(main_hovers)}")
    fields = {f for _, f in main_hovers[0].tooltips}
    check(any("read_name" in f for f in fields), "read tooltip must name the read")
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
        for key in ("lines", "hit", "bin_off", "bin_reads", "read_off", "read_lines"):
            check(key in args, f"{name} callback missing arg {key}")


def test_selection_glyphs(tmp):
    rd = make_region_data(tmp, "sel")
    p = _panels(report_layout(rd))
    r = line_renderer(p["main"])
    check(r.selection_glyph is not None, "selected reads need a distinct glyph")
    check(r.nonselection_glyph is not None, "unselected reads must dim")
    # Alpha is a per-row field so a layer hidden by a Tier-3 toggle stays hidden even
    # when selected; only the nonselection state is a hard zero.
    check(r.glyph.line_alpha == "alpha" and r.selection_glyph.line_alpha == "alpha",
          "base and selected glyphs must honour the per-row alpha field")
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
        check("hit0" in args and "hit1" in args,
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
        for key in ("other", "other_g2l", "other_off", "other_bins"):
            check(key in args, f"{name} hover missing {key} for the partner panel")
        check(args["other"] is hit_renderer(p[partner]).data_source,
              f"{name} must target {partner}'s hit layer")
        check(args["hit"] is hit_renderer(p[name]).data_source,
              f"{name} must stripe its own hovered bin")


def test_nonselected_reads_hidden(tmp):
    rd = make_region_data(tmp, "hid")
    r = line_renderer(_panels(report_layout(rd))["main"])
    check(r.nonselection_glyph.line_alpha == 0.0,
          "reads outside the hovered bin must be invisible, not merely dimmed")


def test_coverage_by_read_transpose(tmp):
    """CoverageBins.by_read must be the exact transpose of the bin index."""
    rd = make_region_data(tmp, "tr")
    for cov in (rd.cov1, rd.cov2):
        off, bins = cov.by_read()
        check(len(off) == len(cov.read_names) + 1, "offsets must be n_reads + 1")
        check(int(off[-1]) == len(cov.read_idx), "must hold every (bin, read) pair")
        for r in range(len(cov.read_names)):
            got = set(bins[off[r]:off[r + 1]].tolist())
            expected = {b for b in range(len(cov.bins))
                        if r in cov.read_idx[cov.bin_offsets[b]:cov.bin_offsets[b + 1]]}
            check(got == expected,
                  f"read {cov.read_names[r]}: bins {got} != {expected}")


def test_band_detection_data(tmp):
    """Hovering a stacked band narrows the highlight to that read class.

    The band is resolved client-side from the cursor's value-axis position against
    this bin's cumulative band tops, so those tops must be on the hit source and the
    read -> band map must agree with the coverage stacking order.
    """
    from qaview.glyphs import _read_bucket_index

    rd = make_region_data(tmp, "band")
    n_series = len(rd.bucket_keys) + 1

    for orient, asm, cov in (("x", rd.asm2, rd.cov2), ("y", rd.asm1, rd.cov1)):
        p = coverage_figure(rd, asm, orient)
        data = hit_renderer(p).data_source.data
        for s in range(n_series):
            check(f"cum{s}" in data, f"{orient}: missing band top cum{s}")
            check(len(data[f"cum{s}"]) == len(cov.bins),
                  f"{orient}: cum{s} must have one value per bin")
        # Bands must be non-decreasing and end at the total.
        for i in range(len(cov.bins)):
            tops = [data[f"cum{s}"][i] for s in range(n_series)]
            check(tops == sorted(tops), f"{orient}: bin {i} band tops not monotonic: {tops}")
        check([round(v, 9) for v in data[f"cum{n_series - 1}"]] ==
              [round(v, 9) for v in cov.total()],
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
    check(args["band_names"] == list(rd.bucket_keys) + ["unassigned"],
          f"band names must mirror the stack order, got {args['band_names']}")
    check(args["n_series"] == n_series, "n_series must match the number of stacked bands")
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
    divs = [c for c in layout.children if isinstance(c, Div)]
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
              [[c, s, e, o] for c, s, e, o in rd.asm1_segs_info],
              "ysegs must match asm1_segs_info exactly")


def main():
    tmp = tempfile.mkdtemp(prefix="qaview_interaction_tests_")
    tests = [test_read_to_line_csr, test_bin_to_line_roundtrip,
             test_coverage_read_names_are_a_subset, test_hit_layer_present_and_invisible,
             test_hover_tools_wired, test_selection_glyphs, test_reverse_link_registered,
             test_draw_order_and_widths, test_hover_links_both_panels,
             test_nonselected_reads_hidden, test_coverage_by_read_transpose,
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
