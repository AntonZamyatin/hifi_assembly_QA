"""Bokeh figure builders for the interactive report.

Read-correspondence lines, marginal coverage and contig bars, wired together with
linked ranges, wheel zoom, a shared crosshair, and hover/tap selection in CustomJS.

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
    ColumnDataSource, CrosshairTool, CustomJS, CustomJSHover, Div, HoverTool, Quad, Range1d,
    Scatter, Segment as SegmentGlyph, Spacer, Span, TapTool, WheelZoomTool,
)
from bokeh.plotting import figure

from qatool.plotters.axis import contig_runs
from qatool.plotters.colors import SUPPORT_COLORS
from qaview.region import LAYER_PARTIAL, MARGIN_BINS, coverage_gaps, marginal_bins

SUPPLEMENTARY_COLOR = "#d62728"
# A one-sided read's partial alignments in the other candidate (region.LAYER_PARTIAL):
# drawn only when a hover or a tap selects the read, in their own colour.
PARTIAL_COLOR = "#ff7f0e"
PARTIAL_ALPHA = 1.0
# A large indel on a partial row: a cross on its line, with it on selection only.
INDEL_COLOR = "#4a148c"
# a read picked in the read browser: rings on its best-placement lines
PICK_COLOR = "#e6007e"
PICK_RING = 18
# coverage dropouts on the marginals; drawn past any coverage, the axis clips them
DROPOUT_COLOR = "#e41a1c"
DROPOUT_TOP = 1e9

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


# Per-layer colour, indexed by qaview.region.LAYER_*. None means "take it from
# SUPPORT_COLORS[class]".
LAYER_COLOR = {
    0: SUPPORT_COLORS["common"],         # common — bottom; matches the marginals
    1: "black",                          # unclassified ("other")
    2: None,                             # support — coloured by class
    3: SUPPLEMENTARY_COLOR,              # supplementary — top
    4: PARTIAL_COLOR,                    # partial: a one-sided read's partial alignments, on selection only
}
LINE_WIDTH = 3.0                         # every read line, in screen pixels: visible at any zoom


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


# ---- line style as arrays ------------------------------------------------------------
# Every per-line column goes to the browser as one binary array: Bokeh encodes a
# Python list element by element, which on a 70,000-line locus was most of the
# view's time and size. Colours are BokehJS's packed uint32, 0xRRGGBBAA, which a
# ColorSpec takes as is. The page restyles lines itself (static/app.js
# `lineStyle`, mirroring `line_style`), so a colouring change needs no server.

LINE_ALPHA = 1.0                         # read lines are opaque: always plainly visible


def pack_rgba(color: str, alpha: int = 255) -> int:
    r, g, b = _hex_rgb(color if color.startswith("#") else to_hex(color))
    return (r << 24) | (g << 16) | (b << 8) | alpha


def ratio_strengths(ratio, C: float) -> np.ndarray:
    """`ratio_strength` over an array; NaN (not scored) gives 0."""
    r = np.asarray(ratio, dtype=float)
    with np.errstate(invalid="ignore"):
        t = np.log1p(np.minimum(np.abs(r), C)) / math.log1p(C)
    t = np.where(np.isnan(r), 0.0, t)
    return np.where(r > 0, t, -t)


def ratio_rgba(t) -> np.ndarray:
    """`ratio_color` over an array, packed."""
    t = np.asarray(t, dtype=float)
    a = np.minimum(np.abs(t), 1.0)[:, None]
    mid = np.array(RATIO_MID, float)
    end = np.where((t > 0)[:, None], np.array(RATIO_POS, float), np.array(RATIO_NEG, float))
    rgb = np.round(mid + (end - mid) * a).astype(np.uint32)
    return (rgb[:, 0] << 24) | (rgb[:, 1] << 16) | (rgb[:, 2] << 8) | np.uint32(255)


def line_style(color_by: str, layer, supp, cls_rgba, ratio, visible_layers=None,
               line_alpha: float = LINE_ALPHA, C: float = 200.0) -> tuple:
    """(colour uint32, alpha float32, selected alpha float32) per line.

    `cls_rgba` is each line's support-class colour, `ratio` its read's log-ratio
    (NaN: not scored). Hidden layers are drawn at alpha 0 rather than removed: the
    hover/selection CSR indexes address rows by position. Partial rows
    (LAYER_PARTIAL) are never drawn unselected, and orange in every mode.
    """
    layer = np.asarray(layer, dtype=np.int64)
    nl = len(LAYER_COLOR)
    base = np.array([pack_rgba(LAYER_COLOR[i]) if LAYER_COLOR[i] else 0 for i in range(nl)],
                    dtype=np.uint32)
    has_base = np.array([LAYER_COLOR[i] is not None for i in range(nl)])
    if color_by == "class":
        color = np.where(has_base[layer], base[layer], np.asarray(cls_rgba, dtype=np.uint32))
    elif color_by == "block":
        # Convention inherited from reads_dotplot.py: black backbone, supplementary red.
        color = np.where(np.asarray(supp, dtype=bool), np.uint32(pack_rgba(SUPPLEMENTARY_COLOR)),
                         np.uint32(pack_rgba("black")))
    elif color_by == "log-ratio":
        ratio = np.asarray(ratio, dtype=float)
        scored = ~np.isnan(ratio)
        t = ratio_strengths(ratio, C)
        color = np.where(scored, ratio_rgba(t), np.uint32(pack_rgba(UNSCORED_COLOR)))
    else:
        raise ValueError("color_by must be 'class', 'block' or 'log-ratio'")
    shown = (np.ones(len(layer), bool) if visible_layers is None
             else np.isin(layer, [int(x) for x in visible_layers]))
    alpha = np.where(shown, line_alpha, 0.0)
    part = layer == 4
    color = np.where(part, np.uint32(pack_rgba(PARTIAL_COLOR)), color)
    sel_alpha = np.where(part, np.where(shown, PARTIAL_ALPHA, 0.0), alpha)
    alpha = np.where(part, 0.0, alpha)
    return color.astype(np.uint32), alpha.astype(np.float32), sel_alpha.astype(np.float32)


def style_spec() -> dict:
    """The constants `line_style` uses, for the page's own restyling."""
    return {"layer_rgba": [pack_rgba(LAYER_COLOR[i]) if LAYER_COLOR[i] else None
                           for i in range(len(LAYER_COLOR))],
            "partial_layer": 4, "partial_rgba": pack_rgba(PARTIAL_COLOR),
            "partial_alpha": PARTIAL_ALPHA,
            "supp_rgba": pack_rgba(SUPPLEMENTARY_COLOR), "black_rgba": pack_rgba("black"),
            "unscored_rgba": pack_rgba(UNSCORED_COLOR), "ratio_mid": list(RATIO_MID),
            "ratio_pos": list(RATIO_POS), "ratio_neg": list(RATIO_NEG),
            "line_alpha": LINE_ALPHA}


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
    """What a coverage marginal stacks, bottom to top: `read_band[g]` is the band
    of global read g (marginal_bins and REBIN_JS bin by it), so whichever stacking
    is drawn, hovering a band isolates exactly the reads counted in it.
    """

    keys: list
    colors: list
    read_band: list


def class_stacking(rd) -> Stacking:
    """Support classes, then the unassigned remainder -- the original marginals."""
    keys = list(rd.bucket_keys) + ["unassigned"]
    colors = [SUPPORT_COLORS[k] for k in rd.bucket_keys] + ["#6e6e6e"]
    return Stacking(keys, colors, _read_bucket_index(rd))


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

    Binned from the same placements as the class stacking, so the bands sum to
    the same total.
    """
    w, s = f"{RATIO_WEAK:g}", f"{RATIO_STRONG:g}"
    keys = [f"indifferent (|r| < {w})",
            f"{rd.asm2} ({w}-{s} nats)", f"{rd.asm2} (>= {s} nats)",
            f"{rd.asm1} ({w}-{s} nats)", f"{rd.asm1} (>= {s} nats)",
            "not scored"]
    colors = [ratio_color(0.0), ratio_color(-0.5), ratio_color(-1.0),
              ratio_color(0.5), ratio_color(1.0), UNSCORED_COLOR]
    read_band = [ratio_band(ratios.get(n)) for n in rd.read_names]
    return Stacking(keys, colors, read_band)


def _read_ratios(rd, ratios) -> np.ndarray:
    """Each read's log-ratio, in `rd.read_names` order; NaN where not scored."""
    if ratios is None:
        return np.full(len(rd.read_names), np.nan)
    return np.array([ratios.get(n, np.nan) for n in rd.read_names], dtype=float)


def _line_sources(rd, color_by: str, visible_layers=None, line_alpha: float = LINE_ALPHA,
                  ratios=None, dropped=None, C: float = 200.0) -> tuple:
    """(lines, reads, tables): the line rows, a per-read table, and the names the
    rows' integer codes stand for.

    `rd.lines` holds single segments (x0,y0)-(x1,y1), so the `segment` glyph fits
    directly and hit-testing works per line -- no MultiLine vertex bookkeeping.
    Rows arrive already sorted back-to-front by `region._order_for_drawing`. A row
    names its read, class and contigs by index (`read` into the `reads` table,
    `cls`, `c1`, `c2` into `tables`); the hover tooltip looks them up.
    """
    n = len(rd.lines)
    read = np.asarray(rd.lines.read_idx, dtype=np.int32)
    layer = np.asarray(rd.lines.layer, dtype=np.int32)
    supp = np.asarray(rd.lines.supplementary, dtype=bool)
    classes = [rd.read_class.get(nm, "other") for nm in rd.read_names]
    cls_names = sorted(set(classes) | {"other"})
    code = {c: i for i, c in enumerate(cls_names)}
    read_cls = np.array([code[c] for c in classes], dtype=np.int32)
    cls_rgba = np.array([pack_rgba(SUPPORT_COLORS.get(c, SUPPORT_COLORS["other"]))
                         for c in cls_names], dtype=np.uint32)
    cls = read_cls[read] if n else np.zeros(0, np.int32)
    ratio = _read_ratios(rd, ratios)
    color, alpha, sel_alpha = line_style(color_by, layer, supp, cls_rgba[cls],
                                         ratio[read] if n else np.zeros(0),
                                         visible_layers, line_alpha, C)
    # what each row is: a full-length placement, or which partial alignment
    labels = list(getattr(rd.lines, "part", None) or [])
    part_names = ["full-length"] + sorted({x for x in labels if x})
    pcode = {x: i for i, x in enumerate(part_names)}
    part = (np.array([pcode[x] if x else 0 for x in labels], dtype=np.int32) if labels
            else np.zeros(n, np.int32))

    def codes(values):
        if not n:
            return [], np.zeros(0, np.int32)
        names, inv = np.unique(np.asarray(values, dtype=object), return_inverse=True)
        return [str(x) for x in names], inv.astype(np.int32)

    c1_names, c1 = codes(rd.lines.contig1)
    c2_names, c2 = codes(rd.lines.contig2)
    i32 = lambda a: np.asarray(a, dtype=np.int32)
    lines = ColumnDataSource(dict(
        x0=np.asarray(rd.lines.x0, float), y0=np.asarray(rd.lines.y0, float),
        x1=np.asarray(rd.lines.x1, float), y1=np.asarray(rd.lines.y1, float),
        read=read, cls=cls, layer=layer, supp=supp.astype(np.int8), c1=c1, c2=c2,
        ref1_0=i32(rd.lines.ref1_0), ref1_1=i32(rd.lines.ref1_1),
        ref2_0=i32(rd.lines.ref2_0), ref2_1=i32(rd.lines.ref2_1),
        pid1=np.asarray(rd.lines.pid1, np.float32), pid2=np.asarray(rd.lines.pid2, np.float32),
        cls_color=cls_rgba[cls], color=color, alpha=alpha, sel_alpha=sel_alpha, part=part,
    ), name="qv-lines")
    reads = ColumnDataSource(dict(
        name=list(rd.read_names), ratio=ratio,
        ratio_txt=[_ratio_text(nm, ratios, dropped) for nm in rd.read_names],
    ), name="qv-reads")
    return lines, reads, {"classes": cls_names, "c1": c1_names, "c2": c2_names, "parts": part_names}


def _line_source(rd, color_by: str, visible_layers=None, line_alpha: float = LINE_ALPHA,
                 ratios=None, dropped=None, C: float = 200.0) -> ColumnDataSource:
    """The line rows alone (see `_line_sources`)."""
    return _line_sources(rd, color_by, visible_layers, line_alpha, ratios, dropped, C)[0]


def line_renderer(fig):
    """The read-line glyph renderer (created first in `dotplot_figure`)."""
    return fig.renderers[0]


def hit_renderer(fig):
    """The invisible per-bin hover target (added last in `coverage_figure`)."""
    return fig.renderers[-1]


# k-mer similarity background (qaview.similarity): one blue hue, light at the
# display floor -- by default the recruitment's pid -- dark at 100%.
BACKGROUND_COLORS = ("#c6dbef", "#08306b")


def _background_index(rd, assembly: str, n: int, bs: int) -> np.ndarray:
    """For each of an axis's n background bins in `rd`'s layout, the bin of the
    built layout (`similarity.axis_map`) its centre comes from -- an arranged axis
    (`region.arrange`) takes the image along, to within a bin."""
    from qaview.similarity import axis_map

    built = {(c, s, e): o for c, s, e, o in axis_map(rd.segments, assembly)}
    segs, flips = rd.segs_info(assembly), rd.flips_of(assembly)
    if not segs:
        return np.arange(n)
    offs = np.array([x[3] for x in segs], float)
    length = np.array([x[2] - x[1] for x in segs], float)
    base = np.array([built.get((c, s, e), o) for c, s, e, o in segs], float)
    centre = (np.arange(n) + 0.5) * bs
    j = np.clip(np.searchsorted(offs, centre, "right") - 1, 0, len(segs) - 1)
    rel = np.clip(centre - offs[j], 0, length[j])
    rel = np.where(np.asarray(flips, bool)[j], length[j] - rel, rel)
    return np.clip(((base[j] + rel) // bs).astype(np.int64), 0, n - 1)


def _draw_background(p, bg: dict, floor: float, rd=None):
    """The cells at or above `floor` as one image under the reads (level underlay):
    one glyph, so drawing it costs the same however many cells light up."""
    from qaview.similarity import image_rgba

    img = image_rgba(bg, floor, BACKGROUND_COLORS)
    if rd is not None and rd.flips:                 # an arranged view: move the cells with it
        bs = bg["binsize"]
        img = img[np.ix_(_background_index(rd, rd.asm1, bg["n1"], bs),
                         _background_index(rd, rd.asm2, bg["n2"], bs))]
    bs = bg["binsize"]
    return p.image_rgba(image=[img], x=0, y=0, dw=bg["n2"] * bs, dh=bg["n1"] * bs,
                        level="underlay", name="qv-background")


def _draw_indels(p, rd, lines):
    """Large indels on one-sided reads' partial rows (RegionData.indels): crosses
    drawn only while their read is selected -- a hover over its marginal bin or a
    tap -- like the rows they sit on. Their text is in the row's tooltip."""
    ind = getattr(rd, "indels", None) or {}
    if not len(ind.get("x", ())):
        return None
    src = ColumnDataSource(dict(x=np.asarray(ind["x"], float), y=np.asarray(ind["y"], float),
                                read=np.asarray(ind["read"], np.int32), text=list(ind["text"])),
                           name="qv-indels")
    mark = dict(x="x", y="y", marker="x", size=12, line_color=INDEL_COLOR, line_width=2.5)
    r = p.scatter(source=src, name="qv-indels-glyph", line_alpha=0.0, **{k: v for k, v in mark.items()
                                                                         if k not in ("x", "y")},
                  x="x", y="y")
    r.selection_glyph = Scatter(**mark, line_alpha=1.0)
    r.nonselection_glyph = Scatter(**mark, line_alpha=0.0)
    lines.selected.js_on_change("indices", CustomJS(args=dict(lines=lines, marks=src), code="""
        const read = lines.data["read"], want = new Set();
        for (const i of lines.selected.indices) want.add(read[i]);
        const mr = marks.data["read"], sel = [];
        for (let j = 0; j < mr.length; j++) if (want.has(mr[j])) sel.push(j);
        marks.selected.indices = sel;
    """))
    return r


def dotplot_figure(rd, color_by: str = "class", line_alpha: float = LINE_ALPHA,
                   width: int = MAIN_WIDTH, height: int | None = None,
                   visible_layers=None, ratios=None, dropped=None, C: float = 200.0,
                   background: dict | None = None, floor: float = 0.997,
                   reads: str = "all"):
    """Read-correspondence dotplot: asm2 on X, asm1 on Y.

    `background` (qaview.similarity) draws the k-mer similarity of the segments
    underneath, from `floor` (the locus's recruitment pid) up.
    `reads="hover"` draws no read until reads are selected -- by hovering a
    marginal bin or tapping -- so a frame costs the background alone; "all"
    draws every read line, as before.
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
    lines, reads_src, tables = _line_sources(rd, color_by, visible_layers, line_alpha,
                                             ratios=ratios, dropped=dropped, C=C)
    r = p.segment(x0="x0", y0="y0", x1="x1", y1="y1", source=lines,
                  color="color", line_width=LINE_WIDTH, alpha="alpha", name="qv-lines-glyph")
    # Selecting reads (by tap, or from a coverage bin) hides every other read outright,
    # so a hovered bin shows only the reads that pass through it.  The selected glyph
    # honours `sel_alpha` -- a row's alpha, so a hidden layer stays hidden even when
    # selected, except partial rows, which only a selection shows.
    r.selection_glyph = SegmentGlyph(x0="x0", y0="y0", x1="x1", y1="y1",
                                     line_color="color", line_width=LINE_WIDTH,
                                     line_alpha="sel_alpha")
    r.nonselection_glyph = SegmentGlyph(x0="x0", y0="y0", x1="x1", y1="y1",
                                        line_color="color", line_width=LINE_WIDTH, line_alpha=0.0)
    if reads == "hover":
        # A scalar alpha of 0 is skipped outright when drawing, not stroked 180,000
        # times; the selection glyph still draws the reads of a hovered bin.
        r.glyph.line_alpha = 0.0
    if background is not None:
        _draw_background(p, background, floor, rd)
    _draw_indels(p, rd, lines)
    # Tapping picks reads: only the read lines answer it, not the hidden indel marks.
    for tool in p.tools:
        if isinstance(tool, TapTool):
            tool.renderers = [r]

    # Rows hold codes; the tooltip turns them back into names (the format names
    # which table).
    names = CustomJSHover(args=dict(reads=reads_src, classes=tables["classes"],
                                    c1=tables["c1"], c2=tables["c2"], parts=tables["parts"]),
                          code="""
        switch (format) {
            case "name": return reads.data["name"][value];
            case "ratio": return reads.data["ratio_txt"][value];
            case "class": return classes[value];
            case "c1": return c1[value];
            case "c2": return c2[value];
            case "yesno": return value ? "yes" : "no";
            case "part": return parts[value];
        }
        return String(value);
        """)
    p.add_tools(HoverTool(
        renderers=[r], line_policy="interp", attachment="vertical",
        tooltips=[
            ("read", "@read{name}"),
            ("class", "@cls{class}"),
            (rd.asm1, "@c1{c1} @ref1_0{0,0}-@ref1_1{0,0}  pid @pid1{0.0000}"),
            (rd.asm2, "@c2{c2} @ref2_0{0,0}-@ref2_1{0,0}  pid @pid2{0.0000}"),
            ("supplementary", "@supp{yesno}"),
            ("log-ratio", "@read{ratio}"),
            ("alignment", "@part{part}"),
        ],
        formatters={"@read": names, "@cls": names, "@c1": names, "@c2": names, "@supp": names,
                    "@part": names}))

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


def _area_columns(left, right, lo, hi) -> dict:
    """A stacked band as one step-shaped area: each bin becomes two points, its
    left and right edge. One area per band rather than one quad per bin: adjacent
    sub-pixel quads leave visible seams between bins."""
    e = np.empty(2 * len(left))
    e[0::2], e[1::2] = left, right
    return dict(e=e, lo=np.repeat(lo, 2), hi=np.repeat(hi, 2))


def _hit_columns(bins: dict, stackings: dict, tops: dict, peak: float) -> dict:
    """The hover target's rows: one per bin, with each band's top as
    `cum_<stacking>_<i>` so a hover can tell which band is under the cursor."""
    total = next(iter(tops.values()))[-1] if tops else np.zeros(len(bins["left"]))
    cols = dict(left=bins["left"], right=bins["right"], seg=bins["seg"], contig=bins["contig"],
                bin_start=bins["bin_start"], bin_end=bins["bin_end"], total=total,
                top=np.full(len(bins["left"]), peak))
    for name, t in tops.items():
        for i, v in enumerate(t):
            cols[f"cum_{name}_{i}"] = v
    return cols


def coverage_figure(rd, assembly: str, orientation: str,
                    width: int = MAIN_WIDTH, height: int | None = None,
                    stacking: "Stacking | None" = None, size: int = MARGIN_SIZE,
                    stackings: dict | None = None, active: str | None = None,
                    value_max: float | None = None):
    """Stacked marginal coverage for one assembly, binned for the whole locus.

    Drawn from `rd.placements`, which `report_layout` re-bins in the page for the
    range in view (REBIN_JS), so zooming in shows finer bins; this draws the first
    frame with the same rule (`region.marginal_bins`).

    `stackings` (name -> Stacking) draws each, only `active` visible, tagged
    `stack:<name>`, so the page can switch between them without the server. Without
    it, `stacking` (default: the support classes) is drawn alone, as "class".
    `value_max` sets the coverage axis's top (default: this marginal's maximum).
    """
    segs = rd.segs_info(assembly)
    extent = rd.plot_w if orientation == "x" else rd.plot_h
    if stackings is None:
        stackings, active = {"class": stacking or class_stacking(rd)}, "class"
    active = active or next(iter(stackings))

    if orientation == "x":
        p = figure(width=width, height=size, x_range=Range1d(0, rd.plot_w),
                   tools="pan,box_zoom,reset")
    else:
        p = figure(width=size, height=height or main_height(rd, width),
                   y_range=Range1d(0, rd.plot_h), tools="pan,box_zoom,reset")

    tops, bins = {}, None
    for name, st in stackings.items():
        bins = marginal_bins(segs, rd.placements[assembly], st.read_band, len(st.keys), 0, extent,
                             flips=rd.flips_of(assembly))
        cum = np.cumsum(bins["series"], axis=0)
        tops[name] = list(cum)
        base = np.zeros(len(bins["left"]))
        for i, (top, color) in enumerate(zip(cum, st.colors)):
            src = ColumnDataSource(_area_columns(bins["left"], bins["right"], base, top))
            kw = dict(source=src, fill_color=color, fill_alpha=1.0, visible=(name == active))
            if orientation == "x":
                r = p.varea(x="e", y1="lo", y2="hi", **kw)
            else:
                r = p.harea(y="e", x1="lo", x2="hi", **kw)
            r.tags = [f"stack:{name}"]          # on the renderer: a glyph method gives it the glyph
            base = top

    total = tops[active][-1] if len(bins["left"]) else np.zeros(0)
    m = value_max if value_max else (float(total.max()) if len(total) else 1.0)
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

    # Coverage dropouts: a thin red mark wherever no read lies, exact to the bp and
    # drawn at least a pixel wide, so a gap narrower than a bin still shows.
    gl, gr = coverage_gaps(segs, rd.placements[assembly])
    gap = dict(fill_color=DROPOUT_COLOR, fill_alpha=0.35, line_color=DROPOUT_COLOR, line_width=1)
    if orientation == "x":
        p.quad(left=gl, right=gr, bottom=0, top=DROPOUT_TOP, **gap)
    else:
        p.quad(left=0, right=DROPOUT_TOP, bottom=gl, top=gr, **gap)

    # Invisible full-height target, one quad per bin.  Hovering the stacked series
    # directly would miss any bucket that is empty at that bin; this always reports the
    # bin under the cursor, and its selection glyph doubles as the bin highlight.
    hit_src = ColumnDataSource(_hit_columns(bins, stackings, tops, m))
    if orientation == "x":
        quad = dict(left="left", right="right", bottom=0, top="top")
        none = dict(left="left", right="right", bottom=0, top=0)
    else:
        quad = dict(left=0, right="top", bottom="left", top="right")
        none = dict(left=0, right=0, bottom="left", top="right")
    hit = p.quad(source=hit_src, fill_alpha=0.0, line_alpha=0.0, fill_color="#ffb300", **quad)
    hit.selection_glyph = Quad(fill_color="#ffb300", fill_alpha=0.30, line_alpha=0.0, **quad)
    hit.nonselection_glyph = Quad(fill_alpha=0.0, line_alpha=0.0, **none)

    p.add_tools(HoverTool(renderers=[hit], attachment="vertical", tooltips=[
        ("bin", "@contig:@bin_start{0,0}-@bin_end{0,0}"),
        ("coverage", "@total{0.0}"),
    ]))
    return p


END_MARK = BAR_SIZE - 3     # contig-end triangle, px
# A segment edge this close to its contig's end counts as reaching it: reads seldom
# align to a contig's last few hundred bases, so recruited segments stop short.
CONTIG_END_TOL = 1000


def contig_end_marks(rd, assembly: str, tol: int = CONTIG_END_TOL) -> list:
    """(plot position, direction, contig, "start" | "end", contig length, bp short)
    for each segment edge within `tol` of an end of its contig.

    Direction is the way into the segment along the axis: +1 at a contig's start
    (the segment runs on to larger positions), -1 at its end -- the other way
    round on a flipped segment (`region.arrange`). Needs `rd.contig_lengths`; a RegionData
    without them (an older pickle) has no marks.
    """
    lengths = (getattr(rd, "contig_lengths", None) or {}).get(assembly, {})
    out = []
    for (contig, start, end, off), f in zip(rd.segs_info(assembly), rd.flips_of(assembly)):
        n = lengths.get(contig)
        if n is None:
            continue
        lo, hi = off, off + end - start            # a flipped segment has its start at hi
        if start <= tol:
            out.append((hi, -1, contig, "start", n, max(0, start)) if f else
                       (lo, 1, contig, "start", n, max(0, start)))
        if end >= n - tol:
            out.append((lo, 1, contig, "end", n, max(0, n - end)) if f else
                       (hi, -1, contig, "end", n, max(0, n - end)))
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
    # The figure's own length, read when it runs: the page may resize the report.
    shift = CustomJS(args=dict(src=src, rng=rng, fig=p, along_x=(orientation == "x"), half=half),
                     code="""
        const px = Math.max(1, along_x ? fig.width : fig.height);
        const per = (rng.end - rng.start) / px;
        const d = src.data;
        d.pos = d.base.map((b, i) => b + d.dir[i] * half * per);
        src.change.emit();
    """)
    rng.js_on_change("start", shift)
    rng.js_on_change("end", shift)
    p.js_on_change("width" if orientation == "x" else "height", shift)
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


def _read_bucket_index(rd) -> list:
    """Global read index -> stacked-band index, matching the coverage series order.

    `coverage_figure` stacks the classes in `rd.bucket_keys` order, then the
    unclassified remainder: band i is class i and the last band is 'unassigned'.
    """
    pos = {k: i for i, k in enumerate(rd.bucket_keys)}
    unassigned = len(rd.bucket_keys)
    return [pos.get(rd.read_class.get(n, ""), unassigned) for n in rd.read_names]


# The marginals' page-side binning. `M` is one marginal: `iv` its placements
# (region.Placements), `segs` its axis map, `hit` its hover target, `areas` its
# band areas per stacking, `first`/`last` the columns of `reads` holding each
# read's placement rows; `reads` also holds `band_<stacking>` per read.
# qvRebin mirrors region.marginal_bins (tests/test_rebin.py runs both).
REBIN_JS = """
function qvNice(span, target) {
    const raw = span / target;
    if (raw <= 1) return 1;
    let mag = 1;
    while (mag * 10 <= raw) mag *= 10;
    for (const f of [1, 2, 5, 10]) if (f * mag >= raw) return f * mag;
}

// Bin M's placements over [a, b]; set its areas and hover rows (all but `top`).
// Returns the peak coverage of what was binned.
function qvRebin(M, reads, a, b, target) {
    const S = M.segs.data, I = M.iv.data;
    const bin = qvNice(b - a, target);
    const left = [], right = [], segk = [], contig = [], gs = [], ge = [];
    const base = new Map();
    for (let k = 0; k < S.off.length; k++) {
        const off = S.off[k], len = S.len[k];
        if (off >= b || off + len <= a) continue;
        const j0 = Math.max(0, Math.floor((a - off) / bin));
        const j1 = Math.min(Math.ceil(len / bin), Math.ceil((b - off) / bin));
        base.set(k, [left.length - j0, j0, j1]);
        for (let j = j0; j < j1; j++) {
            const l = j * bin, r = Math.min(len, l + bin), e = S.start[k] + len;
            left.push(off + l); right.push(off + r); segk.push(k); contig.push(S.contig[k]);
            gs.push(S.flip[k] ? e - r : S.start[k] + l);       // a flipped segment runs e -> s
            ge.push(S.flip[k] ? e - l : S.start[k] + r);
        }
    }
    const nb = left.length, names = Object.keys(M.areas);
    const diffs = names.map((n) => M.areas[n].map(() => new Float64Array(nb + 1)));
    const bandOf = names.map((n) => reads.data["band_" + n]);
    for (let i = 0; i < I.p0.length; i++) {
        const h = base.get(I.seg[i]);
        if (h === undefined) continue;
        const off = S.off[I.seg[i]];
        const lo = Math.max(h[1], Math.floor((I.p0[i] - off) / bin));
        const hi = Math.min(h[2], Math.ceil((I.p1[i] - off) / bin));
        if (lo >= hi) continue;
        const w = I.w[i], g = I.read[i];
        for (let s = 0; s < names.length; s++) {
            const d = diffs[s][bandOf[s][g]];
            d[h[0] + lo] += w;
            d[h[0] + hi] -= w;
        }
    }
    const e = new Float64Array(2 * nb);
    for (let j = 0; j < nb; j++) { e[2 * j] = left[j]; e[2 * j + 1] = right[j]; }
    const twice = (v) => { const o = new Float64Array(2 * nb); for (let j = 0; j < nb; j++) o[2 * j] = o[2 * j + 1] = v[j]; return o; };
    const hit = {left: Float64Array.from(left), right: Float64Array.from(right),
                 seg: Int32Array.from(segk), contig, bin_start: gs, bin_end: ge};
    let total = new Float64Array(nb);
    names.forEach((n, s) => {
        let cum = new Float64Array(nb);
        diffs[s].forEach((d, i) => {
            const top = new Float64Array(nb);
            let run = 0;
            for (let j = 0; j < nb; j++) { run += d[j]; top[j] = cum[j] + (Math.abs(run) < 1e-9 ? 0 : run); }
            M.areas[n][i].data = {e, lo: twice(cum), hi: twice(top)};
            hit["cum_" + n + "_" + i] = top;
            cum = top;
        });
        total = cum;
    });
    hit.total = total;
    M.next = hit;
    let peak = 0;
    for (let j = 0; j < nb; j++) if (total[j] > peak) peak = total[j];
    return peak;
}

// Indices of M's current bins that read set `gs` has a placement in.
function qvBinsOf(M, reads, gs) {
    const L = M.hit.data.left, R = M.hit.data.right, I = M.iv.data;
    const first = reads.data[M.first], last = reads.data[M.last];
    const out = new Set();
    for (const g of gs) {
        for (let i = first[g]; i < last[g]; i++) {
            let lo = 0, hi = R.length;              // first bin ending past p0
            while (lo < hi) { const mid = (lo + hi) >> 1; if (R[mid] > I.p0[i]) hi = mid; else lo = mid + 1; }
            for (let j = lo; j < L.length && L[j] < I.p1[i]; j++) out.add(j);
        }
    }
    return Array.from(out);
}

// The reads with a placement in bin b of M; with band >= 0, only those in it.
function qvReadsInBin(M, b, band, read_band) {
    const h = M.hit.data, I = M.iv.data, gs = new Set();
    for (let i = 0; i < I.p0.length; i++)
        if (I.seg[i] === h.seg[b] && I.p0[i] < h.right[b] && I.p1[i] > h.left[b]
            && (band < 0 || read_band[I.read[i]] === band)) gs.add(I.read[i]);
    return gs;
}

// Read g's line rows. The map is built once per lines data (its `read` column).
function qvRowsOf(lines, g) {
    const read = lines.data["read"];
    const cache = window.qvRowsCache || (window.qvRowsCache = new WeakMap());
    let m = cache.get(read);
    if (!m) {
        m = new Map();
        for (let r = 0; r < read.length; r++) {
            const a = m.get(read[r]);
            if (a) a.push(r); else m.set(read[r], [r]);
        }
        cache.set(read, m);
    }
    return m.get(g) || [];
}

// The reads under the current selected line rows.
function qvSelectedReads(lines) {
    const read = lines.data["read"], gs = new Set();
    for (const r of lines.selected.indices) gs.add(read[r]);
    return gs;
}
"""


def _marginal_model(rd, fig, assembly: str, axis: str) -> dict:
    """The data REBIN_JS needs for one marginal (see its comment)."""
    pl = rd.placements[assembly]
    segs = rd.segs_info(assembly)
    areas = {}
    for r in fig.renderers:
        tag = r.tags[0] if r.tags else ""
        if isinstance(tag, str) and tag.startswith("stack:"):
            areas.setdefault(tag.split(":", 1)[1], []).append(r.data_source)
    return dict(
        iv=ColumnDataSource(dict(p0=pl.p0, p1=pl.p1, w=pl.w, read=pl.read, seg=pl.seg)),
        segs=ColumnDataSource(dict(off=np.array([x[3] for x in segs], float),
                                   len=np.array([x[2] - x[1] for x in segs], float),
                                   contig=[x[0] for x in segs],
                                   start=np.array([x[1] for x in segs], float),
                                   flip=np.array(rd.flips_of(assembly), np.int8))),
        hit=hit_renderer(fig).data_source, areas=areas,
        first=f"first_{axis}", last=f"last_{axis}")


def _reads_source(rd, stackings: dict) -> ColumnDataSource:
    """Per read: its band in each stacking and its placement rows on each axis."""
    p1, p2 = rd.placements[rd.asm1], rd.placements[rd.asm2]
    return ColumnDataSource(dict(
        first_x=p2.first, last_x=p2.last, first_y=p1.first, last_y=p1.last,
        **{f"band_{k}": np.asarray(st.read_band, dtype=np.int32) for k, st in stackings.items()}))


def _wire_rebin(rd, main, shared, lines, margs: list, reads):
    """Re-bin both marginals for the dotplot's ranges on every zoom or pan, at most
    once a frame, and rescale their shared coverage axis to the peak in view."""
    cb = CustomJS(args=dict(xr=main.x_range, yr=main.y_range, shared=shared, lines=lines,
                            mx=margs[0], my=margs[1], reads=reads,
                            wx=float(rd.plot_w), wy=float(rd.plot_h), target=MARGIN_BINS),
                  code=REBIN_JS + """
        if (shared.__qvPending) return;
        shared.__qvPending = true;
        requestAnimationFrame(() => {
            shared.__qvPending = false;
            const px = qvRebin(mx, reads, Math.max(0, xr.start), Math.min(wx, xr.end), target);
            const py = qvRebin(my, reads, Math.max(0, yr.start), Math.min(wy, yr.end), target);
            const peak = Math.max(px, py) || 1;
            for (const M of [mx, my]) {
                M.next.top = new Float64Array(M.next.left.length).fill(peak);
                M.hit.selected.indices = [];
                M.hit.data = M.next;
            }
            shared.bounds = [0, peak];
            shared.start = 0;
            shared.end = peak;
            const gs = qvSelectedReads(lines);       // keep a selected read's bins lit
            if (gs.size) for (const M of [mx, my]) M.hit.selected.indices = qvBinsOf(M, reads, gs);
        });
        """)
    for r in (main.x_range, main.y_range):
        r.js_on_change("start", cb)
        r.js_on_change("end", cb)


def _link_bins_to_reads(cov_fig, marg, other_marg, main, status, orientation,
                        reads, stackings: dict, mode: Div, pick: Div):
    """Hover a coverage bin: stripe that bin, light the matching bins on the other
    assembly, and show only the reads passing through it -- unless the page has
    picked a read (`pick`, see `_pick_read`): that one stays shown until cleared.

    The stack is per band, so the cursor's position along the *value* axis says
    which band is under it -- the hover then narrows to that band alone. Above the
    stack there is no band, and the whole bin is shown. `mode.text` names the
    stacking on show ("class" or "ratio"), which the page may switch.

    The bins change with the zoom, so the reads under one are found by scanning the
    placements (`marg`, see REBIN_JS) -- one pass, well under a frame.
    """
    names = {f"names_{k}": list(st.keys) for k, st in stackings.items()}

    cb = CustomJS(
        args=dict(
            lines=line_renderer(main).data_source, M=marg, O=other_marg, reads=reads,
            status=status, mode=mode, names=names, pick=pick,
            vertical=(orientation == "x"), label=cov_fig.name or "",
        ),
        code=REBIN_JS + """
        if (pick.text) return;
        const hit = M.hit, other = O.hit;
        const hovered = cb_data.index.indices;
        if (!hovered || hovered.length === 0) {
            lines.selected.indices = [];
            hit.selected.indices = [];
            other.selected.indices = [];
            status.text = "";
            return;
        }
        const b = hovered[0];
        const m = mode.text;
        const read_band = reads.data["band_" + m], band_names = names["names_" + m];

        // Which stacked band is the cursor in?  The hit target spans the full value
        // axis, so compare the cursor against this bin's cumulative band tops.
        const geom = cb_data.geometry;
        const v = vertical ? geom.y : geom.x;
        // Empty bands are skipped: at the axis (v <= 0) one would match.
        let band = -1, below = 0;
        for (let s = 0; s < band_names.length; s++) {
            const top = hit.data['cum_' + m + '_' + s][b];
            if (top > below && v <= top) { band = s; break; }
            below = top;
        }

        const gs = qvReadsInBin(M, b, band, read_band);       // narrowed to the band
        const rows = [];
        for (const g of gs) for (const r of qvRowsOf(lines, g)) rows.push(r);
        lines.selected.indices = rows;
        hit.selected.indices = [b];
        other.selected.indices = qvBinsOf(O, reads, gs);

        const where = hit.data['contig'][b] + ':'
                    + hit.data['bin_start'][b].toLocaleString() + '-'
                    + hit.data['bin_end'][b].toLocaleString();
        const cls = band >= 0 ? band_names[band] : 'all bands';
        status.text = label + ' ' + where + ' &nbsp;·&nbsp; <b>' + cls
                    + '</b> &nbsp;·&nbsp; ' + gs.size + ' reads';
        """)
    for tool in cov_fig.tools:
        if isinstance(tool, HoverTool):
            tool.callback = cb


def _link_reads_to_bins(main, margs: list, reads):
    """Selecting read lines (a tap) highlights the coverage bins they contribute to."""
    lines = line_renderer(main).data_source
    cb = CustomJS(args=dict(lines=lines, margs=margs, reads=reads), code=REBIN_JS + """
        const gs = qvSelectedReads(lines);
        for (const M of margs) M.hit.selected.indices = gs.size ? qvBinsOf(M, reads, gs) : [];
    """)
    lines.selected.js_on_change("indices", cb)


def _link_bin_click(cov_fig, marg, reads, mode: Div, stackings: dict, orientation: str):
    """Clicking a marginal bin hands every read in it to the page's read browser, as
    a `qv-bin` DOM event: the bin, the reads (indices into `qv-reads`) and each
    one's band in the stacking on show, with the bands' names and colours."""
    cb = CustomJS(args=dict(M=marg, reads=reads, mode=mode, vertical=(orientation == "x"),
                            label=cov_fig.name or "",
                            names={k: list(st.keys) for k, st in stackings.items()},
                            colors={k: list(st.colors) for k, st in stackings.items()}),
                  code=REBIN_JS + """
        const h = M.hit.data, p = vertical ? cb_obj.x : cb_obj.y;
        let lo = 0, hi = h.right.length;                // first bin ending past p
        while (lo < hi) { const mid = (lo + hi) >> 1; if (h.right[mid] > p) hi = mid; else lo = mid + 1; }
        if (lo >= h.left.length || h.left[lo] > p) return;
        const b = lo, gs = qvReadsInBin(M, b, -1, null);
        const m = mode.text, band = reads.data["band_" + m];
        const list = Array.from(gs);
        document.dispatchEvent(new CustomEvent("qv-bin", {detail: {
            assembly: label, contig: h.contig[b], start: h.bin_start[b], end: h.bin_end[b],
            reads: list, bands: list.map((g) => band[g]),
            band_names: names[m], band_colors: colors[m]}}));
        """)
    cov_fig.js_on_event("tap", cb)


def _pick_read(main, margs: list, reads) -> Div:
    """A hidden model the page sets to "<read index>#<n>" to select that read: its
    lines on the dotplot and its bins on both marginals (it may have no lines --
    a read placed in one candidate only). Empty text clears the selection."""
    pick = Div(text="", visible=False, name="qv-pick")
    lines = line_renderer(main).data_source
    # A read's line is a few pixels long on a large locus: a ring at the middle of
    # each of its best-placement lines shows where it is at any zoom.
    marks = ColumnDataSource(dict(x=[], y=[]), name="qv-pick-marks")
    main.scatter("x", "y", source=marks, marker="circle", size=PICK_RING, fill_alpha=0.0,
                 line_color=PICK_COLOR, line_width=2.5)
    pick.js_on_change("text", CustomJS(
        args=dict(lines=lines, margs=margs, reads=reads, marks=marks, partial=LAYER_PARTIAL),
        code=REBIN_JS + """
        const t = cb_obj.text;
        if (!t) {
            lines.selected.indices = [];
            for (const M of margs) M.hit.selected.indices = [];
            marks.data = {x: [], y: []};
            return;
        }
        const g = parseInt(t), d = lines.data;
        const rows = [], x = [], y = [];
        for (const r of qvRowsOf(lines, g)) {
            rows.push(r);
            if (d.layer[r] === partial) continue;      // rings mark the best placements only
            x.push((d.x0[r] + d.x1[r]) / 2);
            y.push((d.y0[r] + d.y1[r]) / 2);
        }
        lines.selected.indices = rows;
        marks.data = {x, y};
        for (const M of margs) M.hit.selected.indices = qvBinsOf(M, reads, [g]);
        """))
    return pick


def _segs_js(rd) -> dict:
    """Each axis's segments as [contig, start, end, plot offset, flipped], for CustomJS."""
    return dict(ysegs=[[c, s, e, o, f] for (c, s, e, o), f in
                       zip(rd.asm1_segs_info, rd.flips_of(rd.asm1))],
                xsegs=[[c, s, e, o, f] for (c, s, e, o), f in
                       zip(rd.asm2_segs_info, rd.flips_of(rd.asm2))],
                a1=rd.asm1, a2=rd.asm2)


def _mark_alignment(rd, main) -> Div:
    """A hidden model the page sets to a JSON list of alignments {assembly, contig,
    start, end} to mark them: each a translucent band across the dotplot where it
    lies inside the segments (horizontal for asm1, vertical for asm2), else the edge
    of the nearest segment of its contig, on the side it lies past. An alignment on
    a contig the locus does not hold marks nothing."""
    mark = Div(text="", visible=False, name="qv-mark")
    band = ColumnDataSource(dict(left=[], right=[], bottom=[], top=[]), name="qv-mark-band")
    edge = ColumnDataSource(dict(x0=[], y0=[], x1=[], y1=[]), name="qv-mark-edge")
    main.quad(left="left", right="right", bottom="bottom", top="top", source=band,
              fill_color=PICK_COLOR, fill_alpha=0.15, line_color=PICK_COLOR, line_alpha=0.5,
              line_width=1)
    main.segment(x0="x0", y0="y0", x1="x1", y1="y1", source=edge, line_color=PICK_COLOR,
                 line_width=4)
    mark.js_on_change("text", CustomJS(
        args=dict(band=band, edge=edge, W=rd.plot_w, H=rd.plot_h, **_segs_js(rd)), code="""
        const bands = {left: [], right: [], bottom: [], top: []}, edges = {x0: [], y0: [], x1: [], y1: []};
        for (const m of (cb_obj.text ? JSON.parse(cb_obj.text) : [])) {
            const onY = m.assembly === a1;
            if (!onY && m.assembly !== a2) continue;
            const spans = [];
            let near = null, dist = Infinity;
            for (const [c, s, e, o, f] of (onY ? ysegs : xsegs)) {
                if (c !== m.contig) continue;
                const lo = Math.max(s, m.start), hi = Math.min(e, m.end);
                if (lo < hi) { spans.push(f ? [o + e - hi, o + e - lo] : [o + lo - s, o + hi - s]); continue; }
                const before = m.end <= s, d = before ? s - m.end : m.start - e;
                // the segment's start edge is at o, unless the segment is shown flipped
                if (d < dist) { dist = d; near = before !== !!f ? o : o + (e - s); }
            }
            for (const [a, b] of spans) {
                bands.left.push(onY ? 0 : a); bands.right.push(onY ? W : b);
                bands.bottom.push(onY ? a : 0); bands.top.push(onY ? b : H);
            }
            if (!spans.length && near !== null) {
                edges.x0.push(onY ? 0 : near); edges.x1.push(onY ? W : near);
                edges.y0.push(onY ? near : 0); edges.y1.push(onY ? near : H);
            }
        }
        band.data = bands;
        edge.data = edges;
        """))
    return mark


def _coordinate_readout(rd, main):
    """Live contig:position readout for the crosshair, both axes."""
    div = Div(text="<i>hover the dotplot for coordinates</i>",
              styles={"font-family": "monospace", "font-size": "12px",
                      "padding": "2px 0 4px 0"})
    cb = CustomJS(
        args=dict(div=div, **_segs_js(rd)),
        code="""
        function loc(segs, p) {
            for (const [c, s, e, o, f] of segs)
                if (p >= o && p < o + (e - s))
                    return c + ':' + Math.round(f ? e - 1 - (p - o) : s + (p - o)).toLocaleString()
                         + (f ? ' (rev)' : '');
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
                  margin_size: int = MARGIN_SIZE, background: dict | None = None,
                  floor: float = 0.997, reads: str = "all"):
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
                          ratios=ratios, dropped=dropped, C=C, background=background,
                          floor=floor, reads=reads)
    # The marginals follow the reads: log-ratio bands once a score exists, support
    # classes otherwise (including before scoring finishes, where every band would
    # otherwise read "not scored").
    # Both stackings are drawn, one shown, so the page can switch with the colouring.
    stackings = {"class": class_stacking(rd)}
    if ratios is not None:
        stackings["ratio"] = ratio_stacking(rd, ratios)
    active = "ratio" if color_by == "log-ratio" and ratios is not None else "class"
    mode = Div(text=active, visible=False, name="qv-mode")
    # One coverage scale for both marginals, as in report_plot: equal coverage is
    # equally long on either, so the two candidates compare at a glance. The scale
    # follows the zoom: the peak of what is in view (`_wire_rebin`).
    cov_top = coverage_figure(rd, rd.asm2, "x", width=w, size=margin_size,
                              stackings=stackings, active=active)
    cov_right = coverage_figure(rd, rd.asm1, "y", width=w, height=h, size=margin_size,
                                stackings=stackings, active=active)
    peak = max([float(np.max(hit_renderer(f).data_source.data["total"], initial=0.0))
                for f in (cov_top, cov_right)]) or 1.0
    for f in (cov_top, cov_right):
        src = hit_renderer(f).data_source
        src.data["top"] = np.full(len(src.data["left"]), peak)
    shared = Range1d(0, peak, bounds=(0, peak))
    cov_top.y_range = shared
    cov_right.x_range = shared
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

    # Wheel zoom on the dotplot zooms both axes. On a marginal it zooms that
    # marginal's genome axis alone -- the top one X, the right one Y, ranges the
    # dotplot shares -- never its coverage axis, which both marginals share as one
    # scale. maintain_focus=False: once one edge reaches its bound, zooming out
    # keeps widening the other side instead of stopping, so the full locus is always
    # reachable again. With the default the view stuck wherever the first edge hit.
    wheel = WheelZoomTool(maintain_focus=False)
    main.add_tools(wheel)
    main.toolbar.active_scroll = wheel
    for fig, dim in ((cov_top, "width"), (cov_right, "height")):
        axis_wheel = WheelZoomTool(dimensions=dim, maintain_focus=False)
        fig.add_tools(axis_wheel)
        fig.toolbar.active_scroll = axis_wheel

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
    reads = _reads_source(rd, stackings)
    mx = _marginal_model(rd, cov_top, rd.asm2, "x")
    my = _marginal_model(rd, cov_right, rd.asm1, "y")
    pick = _pick_read(main, [mx, my], reads)
    _link_bins_to_reads(cov_top, mx, my, main, status, "x", reads, stackings, mode, pick)
    _link_bins_to_reads(cov_right, my, mx, main, status, "y", reads, stackings, mode, pick)
    _link_reads_to_bins(main, [mx, my], reads)
    _wire_rebin(rd, main, shared, line_renderer(main).data_source, [mx, my], reads)
    _link_bin_click(cov_top, mx, reads, mode, stackings, "x")
    _link_bin_click(cov_right, my, reads, mode, stackings, "y")
    mark = _mark_alignment(rd, main)
    readout = _coordinate_readout(rd, main)

    corner = Spacer(width=BAR_SIZE, height=margin_size, tags=["qv:corner"])
    grid = gridplot(
        [[corner, cov_top, None],
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
    # Tags the page finds the panels by, to resize each in place by its borders
    # (static/app.js `applySizes`); the divs carry what lies beside the dotplot.
    main.tags = ["qv:main"]
    cov_top.tags = ["qv:top"]
    cov_right.tags = ["qv:right"]
    xbar.tags = ["qv:xbar"]
    ybar.tags = ["qv:ybar"]
    for div in (readout, status):
        div.tags = ["qv:line", {"bar": BAR_SIZE}]
        div.width = BAR_SIZE + w + margin_size
        div.styles = {**div.styles, "white-space": "nowrap", "overflow": "hidden",
                      "text-overflow": "ellipsis"}
    return column(readout, status, grid, pick, mark)   # pick, mark: hidden, for the page


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
