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

import math
from dataclasses import dataclass

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
MIN_MAIN = 40               # fit_main's floor, per side
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


def _hex_rgb(h: str) -> tuple:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


# Diverging log-ratio scale on the palette's own hues: asm1 in the red of
# `asm1_only`, asm2 in the green of `asm2_only`, through a neutral grey.
RATIO_POS = _hex_rgb(SUPPORT_COLORS["asm1_only"])
RATIO_NEG = _hex_rgb(SUPPORT_COLORS["asm2_only"])
RATIO_MID = (0xd2, 0xd2, 0xca)
UNSCORED_COLOR = "#e6e6e1"


def ratio_strength(ratio: float, C: float) -> float:
    """Map a read's log-ratio into [-1, 1] on a symmetric log scale, saturating at C.

    Linear would wash the verdict out: at the HG002 seed most reads sit within
    +/-1 nat while the reads that decide it reach -27, so on a 0..200 linear scale
    the decisive ones would look almost neutral. log1p keeps both legible.
    A read that fits only one candidate (+/-inf) saturates.
    """
    if ratio != ratio:                                   # NaN
        return 0.0
    mag = min(abs(ratio), C)
    t = math.log1p(mag) / math.log1p(C)
    return t if ratio > 0 else -t


def ratio_color(t: float) -> str:
    end = RATIO_POS if t > 0 else RATIO_NEG
    a = min(abs(t), 1.0)
    rgb = (round(m + (e - m) * a) for m, e in zip(RATIO_MID, end))
    return "#%02x%02x%02x" % tuple(rgb)


def _ratio_text(name, ratios, dropped) -> str:
    if ratios is None:
        return "not scored yet"
    if name in ratios:
        r = ratios[name]
        if r == math.inf:
            return "+inf (placeable only in asm1)"
        if r == -math.inf:
            return "-inf (placeable only in asm2)"
        return f"{r:+.2f} nats"
    if dropped and name in dropped:
        return "dropped: crosses a segment edge"
    return "not scored"


# Log-ratio bands for the marginals, in nats. Fixed rather than scaled by C: the
# clamp decides how much a read may say, not which way it leans. 10 nats is about
# one substitution at Q40, the natural unit for "this read decides something".
RATIO_WEAK, RATIO_STRONG = 1.0, 10.0


@dataclass(frozen=True)
class Stacking:
    """What a coverage marginal stacks, bottom to top.

    `series[asm]` has one row per band and one column per bin; `read_band[g]` is
    the band of global read g. The hover callback uses both, so whichever stacking
    is drawn, hovering a band isolates exactly the reads counted in it.
    """

    keys: list
    colors: list
    series: dict
    read_band: list


def class_stacking(rd) -> Stacking:
    """Support classes, then the unassigned remainder -- the original marginals."""
    keys = list(rd.bucket_keys) + ["unassigned"]
    colors = [SUPPORT_COLORS[k] for k in rd.bucket_keys] + ["#6e6e6e"]
    series = {rd.asm1: list(rd.cov1.bucket_bins or []) + [rd.cov1.bins],
              rd.asm2: list(rd.cov2.bucket_bins or []) + [rd.cov2.bins]}
    return Stacking(keys, colors, series, _read_bucket_index(rd))


def ratio_band(ratio) -> int:
    """0 indifferent, 1-2 asm2 weak/strong, 3-4 asm1 weak/strong, 5 unscored."""
    if ratio is None or ratio != ratio:
        return 5
    mag = abs(ratio)
    if mag < RATIO_WEAK:
        return 0
    strong = mag >= RATIO_STRONG
    return (2 if strong else 1) if ratio < 0 else (4 if strong else 3)


def ratio_stacking(rd, ratios: dict) -> Stacking:
    """Coverage split by which candidate each read favours, and how strongly.

    Recomputed from the stored per-bin read weights, so it is exact -- the bands
    sum to the same total as the class stacking -- and needs no BAM access.
    """
    w, s = f"{RATIO_WEAK:g}", f"{RATIO_STRONG:g}"
    keys = [f"indifferent (|r| < {w})",
            f"{rd.asm2} ({w}-{s} nats)", f"{rd.asm2} (>= {s} nats)",
            f"{rd.asm1} ({w}-{s} nats)", f"{rd.asm1} (>= {s} nats)",
            "not scored"]
    colors = [ratio_color(0.0), ratio_color(-0.5), ratio_color(-1.0),
              ratio_color(0.5), ratio_color(1.0), UNSCORED_COLOR]
    series = {}
    for asm, cov in ((rd.asm1, rd.cov1), (rd.asm2, rd.cov2)):
        bands = [ratio_band(ratios.get(n)) for n in cov.read_names]
        series[asm] = list(cov.rebucket(bands, len(keys)))
    read_band = [ratio_band(ratios.get(n)) for n in rd.read_names]
    return Stacking(keys, colors, series, read_band)


def _line_source(rd, color_by: str, visible_layers=None,
                 line_alpha: float = 0.55, ratios=None, dropped=None,
                 C: float = 200.0) -> ColumnDataSource:
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
    elif color_by == "log-ratio":
        # Rows cannot be re-sorted to put decisive reads on top: the hover and
        # selection callbacks index rows by position. Width and opacity carry the
        # emphasis instead, so indifferent reads recede and decisive ones stand out.
        colors, widths = [], []
        for n in names:
            if ratios is None or n not in ratios:
                colors.append(UNSCORED_COLOR)
                widths.append(1.0)
            else:
                t = ratio_strength(ratios[n], C)
                colors.append(ratio_color(t))
                widths.append(1.0 + 2.5 * abs(t))
    else:
        raise ValueError("color_by must be 'class', 'block' or 'log-ratio'")
    # Hidden layers are drawn at alpha 0 rather than removed: the hover/selection CSR
    # indexes address rows by position, so dropping rows would silently misalign them.
    alphas = [line_alpha if (visible_layers is None or int(lay) in visible_layers) else 0.0
              for lay in layers]
    if color_by == "log-ratio":
        alphas = [a * (0.3 + 0.7 * abs(ratio_strength(ratios[n], C)))
                  if (ratios is not None and n in ratios) else a * 0.3
                  for a, n in zip(alphas, names)]

    return ColumnDataSource(dict(
        x0=rd.lines.x0, y0=rd.lines.y0, x1=rd.lines.x1, y1=rd.lines.y1,
        read_name=names, read_class=classes, color=colors, width=widths, alpha=alphas,
        supplementary=["yes" if s else "no" for s in rd.lines.supplementary],
        contig1=rd.lines.contig1, contig2=rd.lines.contig2,
        ref1_0=rd.lines.ref1_0, ref1_1=rd.lines.ref1_1,
        ref2_0=rd.lines.ref2_0, ref2_1=rd.lines.ref2_1,
        pid1=rd.lines.pid1, pid2=rd.lines.pid2,
        ratio_txt=[_ratio_text(n, ratios, dropped) for n in names],
    ))


def line_renderer(fig):
    """The read-line glyph renderer (created first in `dotplot_figure`)."""
    return fig.renderers[0]


def hit_renderer(fig):
    """The invisible per-bin hover target (added last in `coverage_figure`)."""
    return fig.renderers[-1]


def dotplot_figure(rd, color_by: str = "class", line_alpha: float = 0.55,
                   width: int = MAIN_WIDTH, height: int | None = None,
                   visible_layers=None, ratios=None, dropped=None, C: float = 200.0):
    """Read-correspondence dotplot: asm2 on X, asm1 on Y.

    Depth is conveyed by alpha overplotting rather than merged line width — merging
    is resolution-dependent and would have to be redone on every zoom.
    """
    # bounds: zoom and pan stop at the locus edges -- there is nothing outside it.
    # The marginals and contig bars share these ranges, so they are bounded too.
    p = figure(width=width, height=height or main_height(rd, width),
               x_range=Range1d(0, rd.plot_w, bounds=(0, rd.plot_w)),
               y_range=Range1d(0, rd.plot_h, bounds=(0, rd.plot_h)),
               tools="pan,box_zoom,reset,save,tap", active_drag="pan",
               background_fill_color="white")
    # Bokeh's level-of-detail mode draws every 10th glyph while zooming or panning
    # once there are more than 2000, so reads vanished mid-zoom and came back only
    # after 0.5 s of stillness. Every line is always drawn instead.
    p.lod_threshold = None
    r = p.segment(x0="x0", y0="y0", x1="x1", y1="y1",
                  source=_line_source(rd, color_by, visible_layers, line_alpha,
                                      ratios=ratios, dropped=dropped, C=C),
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
            ("log-ratio", "@ratio_txt"),
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
                    width: int = MAIN_WIDTH, height: int | None = None,
                    stacking: "Stacking | None" = None, size: int = MARGIN_SIZE):
    """Stacked marginal coverage for one assembly.

    Drawn as explicit quads rather than a stacked area: exact-bp bins have a short
    final bin per segment, which a uniform-width area glyph would misplace.
    """
    cov = rd.cov1 if assembly == rd.asm1 else rd.cov2
    left, right = rd.bin_edges(assembly)
    stacking = stacking or class_stacking(rd)
    series, colors = stacking.series[assembly], stacking.colors

    if orientation == "x":
        p = figure(width=width, height=size, x_range=Range1d(0, rd.plot_w),
                   tools="pan,box_zoom,reset")
    else:
        p = figure(width=size, height=height or main_height(rd, width),
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
    # Ticks by the marginal's size: a small one would otherwise run its labels together.
    n_ticks = max(2, size // 45)
    if orientation == "x":
        p.y_range = Range1d(0, m, bounds=(0, m))
        p.xaxis.visible = False
        p.yaxis.ticker.desired_num_ticks = n_ticks
    else:
        p.x_range = Range1d(0, m, bounds=(0, m))
        p.yaxis.visible = False
        p.xaxis.ticker.desired_num_ticks = n_ticks
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


END_MARK = BAR_SIZE - 3     # contig-end triangle, px
# A segment edge this close to its contig's end counts as reaching it: reads seldom
# align to a contig's last few hundred bases, so recruited segments stop short.
CONTIG_END_TOL = 1000


def contig_end_marks(rd, assembly: str, tol: int = CONTIG_END_TOL) -> list:
    """(plot position, direction, contig, "start" | "end", contig length, bp short)
    for each segment edge within `tol` of an end of its contig.

    Direction is the way into the segment: +1 at a contig's start (the segment runs
    on to larger positions), -1 at its end. Needs `rd.contig_lengths`; a RegionData
    without them (an older pickle) has no marks.
    """
    lengths = (getattr(rd, "contig_lengths", None) or {}).get(assembly, {})
    out = []
    for contig, start, end, off in rd.segs_info(assembly):
        n = lengths.get(contig)
        if n is None:
            continue
        if start <= tol:
            out.append((off, 1, contig, "start", n, max(0, start)))
        if end >= n - tol:
            out.append((off + end - start, -1, contig, "end", n, max(0, n - end)))
    return out


def _draw_end_marks(p, rd, assembly: str, orientation: str, length_px: int):
    """Black triangles at contig ends on a contig bar, pointing into the segment.

    A marker is centred on its point, so each sits half its size inside its
    segment; that is a distance in pixels, so the data offset is recomputed as the
    (linked) range zooms.
    """
    marks = contig_end_marks(rd, assembly)
    if not marks:
        return None
    rng = p.x_range if orientation == "x" else p.y_range
    half = END_MARK / 2
    per_px = (rng.end - rng.start) / max(1, length_px)
    base = [m[0] for m in marks]
    dirs = [m[1] for m in marks]
    # Markers are drawn pointing up; angles turn them counter-clockwise on screen.
    if orientation == "x":
        angle = [-math.pi / 2 if d > 0 else math.pi / 2 for d in dirs]
    else:
        angle = [0.0 if d > 0 else math.pi for d in dirs]
    src = ColumnDataSource(dict(
        pos=[b + d * half * per_px for b, d in zip(base, dirs)], base=base, dir=dirs,
        mid=[0.5] * len(marks), angle=angle,
        label=[f"{c} {which} ({n:,} bp)" + (f", {short:,} bp beyond the segment" if short else "")
               for _, _, c, which, n, short in marks]))
    xy = dict(x="pos", y="mid") if orientation == "x" else dict(x="mid", y="pos")
    r = p.scatter(**xy, source=src, marker="triangle", size=END_MARK, angle="angle",
                  fill_color="#000000", line_color="#ffffff", line_width=1)
    shift = CustomJS(args=dict(src=src, rng=rng, px=max(1, length_px), half=half), code="""
        const per = (rng.end - rng.start) / px;
        const d = src.data;
        d.pos = d.base.map((b, i) => b + d.dir[i] * half * per);
        src.change.emit();
    """)
    rng.js_on_change("start", shift)
    rng.js_on_change("end", shift)
    p.add_tools(HoverTool(renderers=[r], tooltips=[("contig", "@label")]))
    return r


def contig_bar_figure(rd, assembly: str, orientation: str, colors: dict,
                      width: int = MAIN_WIDTH, height: int | None = None, rng=None):
    """Thin coloured bar marking which contig each stretch of the axis belongs to.

    A black triangle at a segment's edge marks the start or end of its contig,
    pointing into the segment. `rng` is the axis range to share (the dotplot's),
    so the triangles follow its zoom.
    """
    segs_info = rd.segs_info(assembly)
    total = rd.plot_w if orientation == "x" else rd.plot_h
    rng = rng if rng is not None else Range1d(0, total)

    if orientation == "x":
        p = figure(width=width, height=BAR_SIZE, x_range=rng,
                   y_range=Range1d(0, 1), tools="")
    else:
        p = figure(width=BAR_SIZE, height=height or main_height(rd, width),
                   y_range=rng, x_range=Range1d(0, 1), tools="")

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
    _draw_end_marks(p, rd, assembly, orientation, p.width if orientation == "x" else p.height)
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


def _link_bins_to_reads(rd, cov_fig, cov, other_fig, other_cov, main, status, orientation,
                        stacking: "Stacking | None" = None):
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
    stacking = stacking or class_stacking(rd)
    n_series = len(stacking.keys)

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
            read_band=stacking.read_band, n_series=n_series,
            band_names=list(stacking.keys),
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


def fit_main(rd, max_width: int = MAIN_WIDTH, max_height: int | None = None) -> tuple:
    """(width, height) of the dotplot: square pixels, inside the given box.

    Width is tried first; if the resulting height overflows `max_height` the plot is
    shrunk to that height instead, so the whole report fits a browser window. The
    floor is small on purpose: a larger one would push the report out of a small
    window rather than shrink it.
    """
    ratio = rd.plot_h / max(rd.plot_w, 1)
    bw = max(MIN_MAIN, int(max_width))
    w = bw
    h = main_height(rd, w)
    if max_height is not None and h > max_height:
        h = max(MIN_MAIN, int(max_height))
        # Never wider than the box: main_height's clamps can make h/ratio exceed it.
        w = min(bw, max(MIN_MAIN, int(round(h / max(ratio, 1e-9)))))
    return w, h


def report_layout(rd, color_by: str = "class", visible_layers=None,
                  ratios=None, dropped=None, C: float = 200.0,
                  max_width: int = MAIN_WIDTH, max_height: int | None = None,
                  margin_size: int = MARGIN_SIZE):
    """Assemble the full interactive report.

    Layout mirrors report_plot: coverage above (asm2) and right (asm1), contig bars
    hugging the dotplot on the bottom (asm2) and left (asm1). `max_width` and
    `max_height` bound the dotplot alone; the marginals add `margin_size` + BAR_SIZE.
    """
    colors = contig_colors(rd)
    # Until a score arrives, log-ratio mode draws the support classes rather than an
    # all-grey plot, matching what the marginals show in the same state.
    line_mode = "class" if color_by == "log-ratio" and ratios is None else color_by
    w, h = fit_main(rd, max_width, max_height)
    main = dotplot_figure(rd, color_by=line_mode, width=w, height=h,
                          visible_layers=visible_layers,
                          ratios=ratios, dropped=dropped, C=C)
    # The marginals follow the reads: log-ratio bands once a score exists, support
    # classes otherwise (including before scoring finishes, where every band would
    # otherwise read "not scored").
    stacking = (ratio_stacking(rd, ratios) if color_by == "log-ratio" and ratios is not None
                else class_stacking(rd))
    cov_top = coverage_figure(rd, rd.asm2, "x", width=w, stacking=stacking,
                              size=margin_size)
    cov_right = coverage_figure(rd, rd.asm1, "y", width=w, height=h, stacking=stacking,
                                size=margin_size)
    xbar = contig_bar_figure(rd, rd.asm2, "x", colors, width=w, rng=main.x_range)
    ybar = contig_bar_figure(rd, rd.asm1, "y", colors, width=w, height=h, rng=main.y_range)

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
    # The read tooltip starts switched off: it covered the shift+drag rectangle.
    # It stays in the dotplot's toolbar, one click away.
    main.toolbar.active_inspect = [crosshair]

    # Wheel zoom on the dotplot only, so scrolling a marginal cannot desync the axes.
    # maintain_focus=False: once one edge reaches its bound, zooming out keeps
    # widening the other side instead of stopping, so the full locus is always
    # reachable again. With the default the view stuck wherever the first edge hit.
    wheel = WheelZoomTool(maintain_focus=False)
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
    _link_bins_to_reads(rd, cov_top, rd.cov2, cov_right, rd.cov1, main, status, "x", stacking)
    _link_bins_to_reads(rd, cov_right, rd.cov1, cov_top, rd.cov2, main, status, "y", stacking)
    _link_reads_to_bins(rd, main, [(cov_top, rd.cov2), (cov_right, rd.cov1)])
    readout = _coordinate_readout(rd, main)

    grid = gridplot(
        [[Spacer(width=BAR_SIZE, height=margin_size), cov_top, None],
         [ybar,                                       main,    cov_right],
         [None,                                       xbar,    None]],
        toolbar_location=None, merge_tools=False,
    )
    # The dotplot's toolbar is shown at the grid's right edge, past the right
    # marginal, so the marginal sits flush against the dotplot. It is the same
    # Toolbar model, so its tools still act on the dotplot alone.
    main.toolbar_location = None
    grid.toolbar = main.toolbar
    grid.toolbar_location = "right"
    # One line each, as wide as the plots: a wrapped line would make the report
    # taller than the cell it was sized for (server.PAD_H counts one line each).
    for div in (readout, status):
        div.width = BAR_SIZE + w + margin_size
        div.styles = {**div.styles, "white-space": "nowrap", "overflow": "hidden",
                      "text-overflow": "ellipsis"}
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
