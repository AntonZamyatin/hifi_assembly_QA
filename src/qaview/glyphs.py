"""Bokeh figure builders for the interactive report.

Phase 1: the static panels — read-correspondence lines, marginal coverage, contig
bars — wired together with linked ranges, wheel zoom and a shared crosshair.
No widgets and no hover/selection yet (Phases 2-3).

Every number comes from `RegionData`; nothing here touches a BAM or a DB.  Contig
colours are taken from the same matplotlib colormap `report_plot` uses so the
interactive view and the published figure stay visually comparable.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

from bokeh.layouts import column, gridplot
from bokeh.models import (
    ColumnDataSource, CrosshairTool, CustomJS, Div, HoverTool, Quad, Range1d,
    Segment as SegmentGlyph, Spacer, Span, TapTool, WheelZoomTool,
)
from bokeh.plotting import figure

from qatool.plotters.axis import contig_runs
from qatool.plotters.colors import SUPPORT_COLORS

SUPPLEMENTARY_COLOR = "#d62728"

# Panel sizes in px.  The matplotlib report uses a 16in figure with 2in marginals and
# 0.15in contig bars; at 100 px/in those become the values below.
MAIN_WIDTH = 1100
MARGIN_SIZE = 200
BAR_SIZE = 15


def main_height(rd, width: int = MAIN_WIDTH) -> int:
    """Dotplot height giving square pixels, as report_plot does with its figsize.

    Without this a region whose two assemblies differ in length renders its diagonals
    off 45 degrees, which is exactly the signal the plot exists to show.  Clamped so a
    very lopsided region stays usable.
    """
    ratio = rd.plot_h / max(rd.plot_w, 1)
    return int(max(200, min(4000, round(width * ratio))))


def contig_colors(rd) -> dict:
    """contig -> hex colour, matching `_draw_contig_bars` in report.py."""
    all_contigs = sorted({c for c, _, _, _ in rd.asm1_segs_info + rd.asm2_segs_info})
    cmap = plt.get_cmap("tab10" if len(all_contigs) <= 10 else "tab20")
    return {c: to_hex(cmap(i % cmap.N)) for i, c in enumerate(all_contigs)}


# Per-layer style, indexed by qaview.region.LAYER_*.  Widths decrease going up so the
# wide grey common backbone stays readable under the narrower, more informative lines.
# A colour of None means "take it from SUPPORT_COLORS[class]".
LAYER_STYLE = {
    0: (SUPPORT_COLORS["common"], 4.0),  # common — bottom, widest; matches the marginals
    1: ("black", 2.0),                   # unclassified ("other")
    2: (None, 1.5),                      # support — coloured by class
    3: (SUPPLEMENTARY_COLOR, 1.5),       # supplementary — top
}


def _line_source(rd, color_by: str, visible_layers=None,
                 line_alpha: float = 0.55) -> ColumnDataSource:
    """One row per read-correspondence line segment.

    `rd.lines` holds single segments (x0,y0)-(x1,y1), so the `segment` glyph fits
    directly and hit-testing works per line — no MultiLine vertex bookkeeping.
    Rows arrive already sorted back-to-front by `region._order_for_drawing`.
    """
    names = [rd.read_names[i] for i in rd.lines.read_idx]
    classes = [rd.read_class.get(n, "other") for n in names]
    layers = list(rd.lines.layer)
    if color_by == "class":
        colors, widths = [], []
        for cls, lay in zip(classes, layers):
            base, w = LAYER_STYLE[int(lay)]
            colors.append(base if base is not None
                          else SUPPORT_COLORS.get(cls, SUPPORT_COLORS["other"]))
            widths.append(w)
    elif color_by == "block":
        # Convention inherited from reads_dotplot.py: black backbone, supplementary red.
        colors = [SUPPLEMENTARY_COLOR if s else "black" for s in rd.lines.supplementary]
        widths = [LAYER_STYLE[int(lay)][1] for lay in layers]
    else:
        raise ValueError("color_by must be 'class' or 'block'")
    # Hidden layers are drawn at alpha 0 rather than removed: the hover/selection CSR
    # indexes address rows by position, so dropping rows would silently misalign them.
    alphas = [line_alpha if (visible_layers is None or int(lay) in visible_layers) else 0.0
              for lay in layers]

    return ColumnDataSource(dict(
        x0=rd.lines.x0, y0=rd.lines.y0, x1=rd.lines.x1, y1=rd.lines.y1,
        read_name=names, read_class=classes, color=colors, width=widths, alpha=alphas,
        supplementary=["yes" if s else "no" for s in rd.lines.supplementary],
        contig1=rd.lines.contig1, contig2=rd.lines.contig2,
        ref1_0=rd.lines.ref1_0, ref1_1=rd.lines.ref1_1,
        ref2_0=rd.lines.ref2_0, ref2_1=rd.lines.ref2_1,
        pid1=rd.lines.pid1, pid2=rd.lines.pid2,
    ))


def line_renderer(fig):
    """The read-line glyph renderer (created first in `dotplot_figure`)."""
    return fig.renderers[0]


def hit_renderer(fig):
    """The invisible per-bin hover target (added last in `coverage_figure`)."""
    return fig.renderers[-1]


def dotplot_figure(rd, color_by: str = "class", line_alpha: float = 0.55,
                   width: int = MAIN_WIDTH, height: int | None = None,
                   visible_layers=None):
    """Read-correspondence dotplot: asm2 on X, asm1 on Y.

    Depth is conveyed by alpha overplotting rather than merged line width — merging
    is resolution-dependent and would have to be redone on every zoom.
    """
    p = figure(width=width, height=height or main_height(rd, width),
               x_range=Range1d(0, rd.plot_w), y_range=Range1d(0, rd.plot_h),
               tools="pan,box_zoom,reset,save,tap", active_drag="pan",
               background_fill_color="white")
    r = p.segment(x0="x0", y0="y0", x1="x1", y1="y1",
                  source=_line_source(rd, color_by, visible_layers, line_alpha),
                  color="color", line_width="width", alpha="alpha")
    # Selecting reads (by tap, or from a coverage bin) hides every other read outright,
    # so a hovered bin shows only the reads that pass through it.  The selected glyph
    # also honours `alpha`, so a hidden layer stays hidden even when it is selected.
    r.selection_glyph = SegmentGlyph(x0="x0", y0="y0", x1="x1", y1="y1",
                                     line_color="color", line_width="width",
                                     line_alpha="alpha")
    r.nonselection_glyph = SegmentGlyph(x0="x0", y0="y0", x1="x1", y1="y1",
                                        line_color="color", line_width="width", line_alpha=0.0)

    p.add_tools(HoverTool(
        renderers=[r], line_policy="interp", attachment="vertical",
        tooltips=[
            ("read", "@read_name"),
            ("class", "@read_class"),
            (rd.asm1, "@contig1 @ref1_0{0,0}-@ref1_1{0,0}  pid @pid1{0.0000}"),
            (rd.asm2, "@contig2 @ref2_0{0,0}-@ref2_1{0,0}  pid @pid2{0.0000}"),
            ("supplementary", "@supplementary"),
        ]))

    # Segment boundaries: pale for same-contig, full strength across contigs.
    for segs_info, dim in ((rd.asm2_segs_info, "height"), (rd.asm1_segs_info, "width")):
        for i in range(1, len(segs_info)):
            same = segs_info[i][0] == segs_info[i - 1][0]
            p.add_layout(Span(location=segs_info[i][3], dimension=dim,
                              line_color="black", line_alpha=0.5 if same else 1.0,
                              line_width=0.5 if same else 1.0))

    p.xaxis.visible = False
    p.yaxis.visible = False
    p.grid.visible = False
    return p


def coverage_figure(rd, assembly: str, orientation: str,
                    width: int = MAIN_WIDTH, height: int | None = None):
    """Stacked marginal coverage for one assembly.

    Drawn as explicit quads rather than a stacked area: exact-bp bins have a short
    final bin per segment, which a uniform-width area glyph would misplace.
    """
    cov = rd.cov1 if assembly == rd.asm1 else rd.cov2
    left, right = rd.bin_edges(assembly)
    series = list(cov.bucket_bins or []) + [cov.bins]
    colors = [SUPPORT_COLORS[k] for k in rd.bucket_keys] + ["#6e6e6e"]

    if orientation == "x":
        p = figure(width=width, height=MARGIN_SIZE, x_range=Range1d(0, rd.plot_w),
                   tools="pan,box_zoom,reset")
    else:
        p = figure(width=MARGIN_SIZE, height=height or main_height(rd, width),
                   y_range=Range1d(0, rd.plot_h), tools="pan,box_zoom,reset")

    # Each bin becomes two points (its left and right edge) so one filled area per
    # bucket renders as a step.  Drawn as a single area rather than one quad per bin:
    # adjacent sub-pixel quads leave visible seams between bins at this bin count.
    edges = np.empty(2 * len(left), dtype=float)
    edges[0::2] = left
    edges[1::2] = right

    def _double(vals):
        out = np.empty(2 * len(vals), dtype=float)
        out[0::2] = vals
        out[1::2] = vals
        return out

    base = [0.0] * len(left)
    cums = []                      # cumulative top of each band, per bin
    for values, color in zip(series, colors):
        top = [b + v for b, v in zip(base, values)]
        lo, hi = _double(base), _double(top)
        if orientation == "x":
            p.varea(x=edges, y1=lo, y2=hi, fill_color=color, fill_alpha=1.0)
        else:
            p.harea(y=edges, x1=lo, x2=hi, fill_color=color, fill_alpha=1.0)
        base = top
        cums.append(list(top))

    m = max(base) if base else 1.0
    if orientation == "x":
        p.y_range = Range1d(0, m)
        p.xaxis.visible = False
    else:
        p.x_range = Range1d(0, m)
        p.yaxis.visible = False
    p.grid.visible = False

    # Invisible full-height target, one quad per bin.  Hovering the stacked series
    # directly would miss any bucket that is empty at that bin; this always reports the
    # bin under the cursor, and its selection glyph doubles as the bin highlight.
    n_reads = [int(cov.bin_offsets[i + 1] - cov.bin_offsets[i]) for i in range(len(left))]
    hit_src = ColumnDataSource(dict(
        left=left, right=right, contig=cov.bin_contig,
        bin_start=cov.bin_start, bin_end=cov.bin_end,
        total=base, n_reads=n_reads, bin_index=list(range(len(left))),
        # Band boundaries, so a hover can resolve which stacked class it is over.
        **{f"cum{i}": c for i, c in enumerate(cums)},
    ))
    if orientation == "x":
        hit = p.quad(left="left", right="right", bottom=0, top=m, source=hit_src,
                     fill_alpha=0.0, line_alpha=0.0, fill_color="#ffb300")
        sel = Quad(left="left", right="right", bottom=0, top=m,
                   fill_color="#ffb300", fill_alpha=0.30, line_alpha=0.0)
    else:
        hit = p.quad(left=0, right=m, bottom="left", top="right", source=hit_src,
                     fill_alpha=0.0, line_alpha=0.0, fill_color="#ffb300")
        sel = Quad(left=0, right=m, bottom="left", top="right",
                   fill_color="#ffb300", fill_alpha=0.30, line_alpha=0.0)
    hit.selection_glyph = sel
    hit.nonselection_glyph = Quad(left="left", right="right", bottom=0, top=0,
                                  fill_alpha=0.0, line_alpha=0.0) if orientation == "x" else \
        Quad(left=0, right=0, bottom="left", top="right", fill_alpha=0.0, line_alpha=0.0)

    p.add_tools(HoverTool(renderers=[hit], attachment="vertical", tooltips=[
        ("bin", "@contig:@bin_start{0,0}-@bin_end{0,0}"),
        ("coverage", "@total{0.0}"),
        ("reads", "@n_reads"),
    ]))
    return p


def contig_bar_figure(rd, assembly: str, orientation: str, colors: dict,
                      width: int = MAIN_WIDTH, height: int | None = None):
    """Thin coloured bar marking which contig each stretch of the axis belongs to."""
    segs_info = rd.segs_info(assembly)
    total = rd.plot_w if orientation == "x" else rd.plot_h

    if orientation == "x":
        p = figure(width=width, height=BAR_SIZE, x_range=Range1d(0, total),
                   y_range=Range1d(0, 1), tools="")
    else:
        p = figure(width=BAR_SIZE, height=height or main_height(rd, width),
                   y_range=Range1d(0, total), x_range=Range1d(0, 1), tools="")

    for contig, off, end in contig_runs(segs_info, total):
        if orientation == "x":
            p.quad(left=[off], right=[end], bottom=[0], top=[1],
                   fill_color=colors[contig], line_color=None)
        else:
            p.quad(left=[0], right=[1], bottom=[off], top=[end],
                   fill_color=colors[contig], line_color=None)

    p.xaxis.visible = False
    p.yaxis.visible = False
    p.grid.visible = False
    p.toolbar_location = None
    p.min_border = 0
    return p


def _global_index(rd, cov):
    """(local -> global, global -> local) read-index maps between a CoverageBins and rd.

    `CoverageBins.read_names` is per-assembly and holds only reads that reached a bin,
    while line rows index `RegionData.read_names`.  Precomputing both directions in
    Python keeps the hover callback to plain array lookups.
    """
    pos = {n: i for i, n in enumerate(rd.read_names)}
    l2g = [pos.get(n, -1) for n in cov.read_names]
    g2l = [-1] * len(rd.read_names)
    for local, g in enumerate(l2g):
        if g >= 0:
            g2l[g] = local
    return l2g, g2l


def _read_bucket_index(rd) -> list:
    """Global read index -> stacked-band index, matching the coverage series order.

    `coverage_figure` stacks `bucket_bins` (in `rd.bucket_keys` order) then the
    unbucketed remainder, so band i is bucket i and the last band is 'unassigned'.
    """
    pos = {k: i for i, k in enumerate(rd.bucket_keys)}
    unassigned = len(rd.bucket_keys)
    return [pos.get(rd.read_class.get(n, ""), unassigned) for n in rd.read_names]


def _link_bins_to_reads(rd, cov_fig, cov, other_fig, other_cov, main, status, orientation):
    """Hover a coverage bin: stripe that bin, light the matching bins on the other
    assembly, and show only the reads passing through it.

    The stack is per-class, so the cursor's position along the *value* axis says which
    class is under it — the hover then narrows to that class alone.  Above the stack
    there is no band, and the whole bin is shown.

    Pure client-side array lookups over precomputed indexes — bin -> reads and
    read -> bins on each assembly, read -> line rows, and read -> band — so it stays
    instant over a remote connection.
    """
    rl_off, rl_idx = rd.lines.by_read(len(rd.read_names))
    l2g, _ = _global_index(rd, cov)
    _, other_g2l = _global_index(rd, other_cov)
    other_off, other_bins = other_cov.by_read()
    n_series = len(rd.bucket_keys) + 1

    cb = CustomJS(
        args=dict(
            lines=line_renderer(main).data_source,
            hit=hit_renderer(cov_fig).data_source,
            other=hit_renderer(other_fig).data_source,
            status=status,
            bin_off=cov.bin_offsets.tolist(), bin_reads=cov.read_idx.tolist(),
            l2g=l2g, other_g2l=other_g2l,
            other_off=other_off.tolist(), other_bins=other_bins.tolist(),
            read_off=rl_off.tolist(), read_lines=rl_idx.tolist(),
            read_band=_read_bucket_index(rd), n_series=n_series,
            band_names=list(rd.bucket_keys) + ["unassigned"],
            vertical=(orientation == "x"), label=cov_fig.name or "",
        ),
        code="""
        const hovered = cb_data.index.indices;
        if (!hovered || hovered.length === 0) {
            lines.selected.indices = [];
            hit.selected.indices = [];
            other.selected.indices = [];
            status.text = "";
            return;
        }
        const b = hovered[0];

        // Which stacked band is the cursor in?  The hit target spans the full value
        // axis, so compare the cursor against this bin's cumulative band tops.
        const geom = cb_data.geometry;
        const v = vertical ? geom.y : geom.x;
        let band = -1;
        for (let s = 0; s < n_series; s++) {
            if (v <= hit.data['cum' + s][b]) { band = s; break; }
        }

        const rows = [];
        const partner = new Set();
        let n = 0;
        for (let k = bin_off[b]; k < bin_off[b + 1]; k++) {
            const g = l2g[bin_reads[k]];
            if (g < 0) continue;
            if (band >= 0 && read_band[g] !== band) continue;   // narrow to the class
            n += 1;
            for (let j = read_off[g]; j < read_off[g + 1]; j++) rows.push(read_lines[j]);
            const ol = other_g2l[g];
            if (ol >= 0)
                for (let j = other_off[ol]; j < other_off[ol + 1]; j++)
                    partner.add(other_bins[j]);
        }
        lines.selected.indices = rows;
        hit.selected.indices = [b];
        other.selected.indices = Array.from(partner);

        const where = hit.data['contig'][b] + ':'
                    + hit.data['bin_start'][b].toLocaleString() + '-'
                    + hit.data['bin_end'][b].toLocaleString();
        const cls = band >= 0 ? band_names[band] : 'all classes';
        status.text = label + ' ' + where + ' &nbsp;·&nbsp; <b>' + cls
                    + '</b> &nbsp;·&nbsp; ' + n + ' reads';
        """)
    for tool in cov_fig.tools:
        if isinstance(tool, HoverTool):
            tool.callback = cb


def _link_reads_to_bins(rd, main, cov_figs):
    """Tapping a read line highlights the coverage bins it contributes to."""
    args = dict(lines=line_renderer(main).data_source)
    code_parts = []
    for i, (cov_fig, cov) in enumerate(cov_figs):
        args[f"hit{i}"] = hit_renderer(cov_fig).data_source
        args[f"boff{i}"] = cov.bin_offsets.tolist()
        args[f"breads{i}"] = cov.read_idx.tolist()
        args[f"cnames{i}"] = cov.read_names
        code_parts.append(f"""
        {{
            const local = new Set();
            for (let i = 0; i < cnames{i}.length; i++)
                if (want.has(cnames{i}[i])) local.add(i);
            const bins = [];
            for (let b = 0; b + 1 < boff{i}.length; b++) {{
                for (let k = boff{i}[b]; k < boff{i}[b + 1]; k++) {{
                    if (local.has(breads{i}[k])) {{ bins.push(b); break; }}
                }}
            }}
            hit{i}.selected.indices = bins;
        }}""")
    args["names"] = list(line_renderer(main).data_source.data["read_name"])
    cb = CustomJS(args=args, code="""
        const rows = lines.selected.indices;
        const want = new Set();
        for (const r of rows) want.add(names[r]);
        if (want.size === 0) { """ + "".join(
            f"hit{i}.selected.indices = []; " for i in range(len(cov_figs))
        ) + """ return; }
    """ + "".join(code_parts))
    line_renderer(main).data_source.selected.js_on_change("indices", cb)


def _coordinate_readout(rd, main):
    """Live contig:position readout for the crosshair, both axes."""
    div = Div(text="<i>hover the dotplot for coordinates</i>",
              styles={"font-family": "monospace", "font-size": "12px",
                      "padding": "2px 0 4px 0"})
    cb = CustomJS(
        args=dict(div=div,
                  ysegs=[[c, s, e, o] for c, s, e, o in rd.asm1_segs_info],
                  xsegs=[[c, s, e, o] for c, s, e, o in rd.asm2_segs_info],
                  a1=rd.asm1, a2=rd.asm2),
        code="""
        function loc(segs, p) {
            for (const [c, s, e, o] of segs)
                if (p >= o && p < o + (e - s))
                    return c + ':' + Math.round(s + (p - o)).toLocaleString();
            return '&mdash;';
        }
        div.text = a1 + ' ' + loc(ysegs, cb_obj.y) + ' &nbsp;|&nbsp; '
                 + a2 + ' ' + loc(xsegs, cb_obj.x);
        """)
    main.js_on_event("mousemove", cb)
    return div


def report_layout(rd, color_by: str = "class", visible_layers=None):
    """Assemble the full interactive report.

    Layout mirrors report_plot: coverage above (asm2) and right (asm1), contig bars
    hugging the dotplot on the bottom (asm2) and left (asm1).
    """
    colors = contig_colors(rd)
    h = main_height(rd, MAIN_WIDTH)
    main = dotplot_figure(rd, color_by=color_by, height=h, visible_layers=visible_layers)
    cov_top = coverage_figure(rd, rd.asm2, "x")
    cov_right = coverage_figure(rd, rd.asm1, "y", height=h)
    xbar = contig_bar_figure(rd, rd.asm2, "x", colors)
    ybar = contig_bar_figure(rd, rd.asm1, "y", colors, height=h)

    # Linked ranges: every panel follows the dotplot's viewport.
    cov_top.x_range = main.x_range
    xbar.x_range = main.x_range
    cov_right.y_range = main.y_range
    ybar.y_range = main.y_range

    # One CrosshairTool instance shared by several figures gives a linked cursor:
    # the same Span overlays are rendered on each, so the ruler projects onto the
    # marginals.  Verified against bokeh 3.9.
    crosshair = CrosshairTool(
        overlay=(Span(dimension="width", line_dash="dashed", line_color="#888", line_width=1),
                 Span(dimension="height", line_dash="dashed", line_color="#888", line_width=1)),
    )
    for p in (main, cov_top, cov_right):
        p.add_tools(crosshair)

    # Wheel zoom on the dotplot only, so scrolling a marginal cannot desync the axes.
    wheel = WheelZoomTool()
    main.add_tools(wheel)
    main.toolbar.active_scroll = wheel

    # Only the dotplot keeps a toolbar.  Merging them instead would make bokeh pick a
    # single active drag/scroll tool across panels with different ones, which both
    # warns and silently changes which panel responds to the wheel.
    for p in (cov_top, cov_right):
        p.toolbar_location = None

    # Interaction wiring (Phase 2).
    cov_top.name = rd.asm2
    cov_right.name = rd.asm1
    status = Div(text="", styles={"font-family": "monospace", "font-size": "12px",
                                  "padding": "0 0 4px 0"})
    _link_bins_to_reads(rd, cov_top, rd.cov2, cov_right, rd.cov1, main, status, "x")
    _link_bins_to_reads(rd, cov_right, rd.cov1, cov_top, rd.cov2, main, status, "y")
    _link_reads_to_bins(rd, main, [(cov_top, rd.cov2), (cov_right, rd.cov1)])
    readout = _coordinate_readout(rd, main)

    grid = gridplot(
        [[Spacer(width=BAR_SIZE, height=MARGIN_SIZE), cov_top, None],
         [ybar,                                       main,    cov_right],
         [None,                                       xbar,    None]],
        toolbar_location=None, merge_tools=False,
    )
    return column(readout, status, grid)


def save_html(layout, path: str, title: str = "qaview", inline: bool = True) -> str:
    """Write a standalone HTML file.

    inline=True embeds bokeh's JS (~1.2 MB, fixed) so the page opens with no network;
    inline=False pulls it from the CDN, giving a ~50 KB file that needs internet.

    Deliberately does not call `output_file`: that sets global bokeh state and would
    silently redirect any later `show()` in the same session to this path.
    """
    from bokeh.plotting import save
    from bokeh.resources import CDN, INLINE
    return save(layout, filename=path, resources=INLINE if inline else CDN, title=title)
