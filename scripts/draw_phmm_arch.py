#!/usr/bin/env python
"""Render the object/call scheme of `qamodel.phmm` to a figure.

Kept as a script rather than drawn once by hand so it can be regenerated when
the API moves. Usage: venv/bin/python scripts/draw_phmm_arch.py out.png
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

SERIF, MONO = "DejaVu Serif", "DejaVu Sans Mono"
INK = "#22252a"
FILL = {"in": "#eef3f8", "api": "#e7f0e6", "eng": "#fdf0e4", "out": "#f2ecf6"}
EDGE = {"in": "#7f9bb5", "api": "#7ba076", "eng": "#d2a273", "out": "#9d84b0"}


def box(ax, x, y, w, h, title, lines, kind, title_mono=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.6,rounding_size=1.2",
                                facecolor=FILL[kind], edgecolor=EDGE[kind],
                                linewidth=1.3, zorder=2))
    ax.text(x + w / 2, y + h - 2.0, title, ha="center", va="top", zorder=3,
            fontname=MONO if title_mono else SERIF, fontsize=10.5,
            fontweight="bold", color=INK)
    for i, ln in enumerate(lines):
        mono = ln.startswith("  ") or "(" in ln or "[" in ln
        ax.text(x + 1.6, y + h - 5.4 - i * 2.7, ln, ha="left", va="top", zorder=3,
                fontname=MONO if mono else SERIF, fontsize=8.2, color="#3c4149")
    return (x + w / 2, y, x + w / 2, y + h)


def arrow(ax, p, q, label=None, rad=0.0, style="-|>", color="#5d646e", lp=0.5):
    ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f"arc3,rad={rad}",
                                 arrowstyle=style, mutation_scale=11,
                                 linewidth=1.1, color=color, zorder=1,
                                 shrinkA=2, shrinkB=3))
    if label:
        mx = p[0] + (q[0] - p[0]) * lp + rad * (q[1] - p[1]) * 0.5
        my = p[1] + (q[1] - p[1]) * lp - rad * (q[0] - p[0]) * 0.5
        ax.text(mx, my, label, ha="center", va="center", fontname=MONO,
                fontsize=7.6, color="#5d646e", zorder=4,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="none", alpha=0.92))


def draw(path):
    fig = plt.figure(figsize=(15.5, 11.2))
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.94])
    ax.set_xlim(0, 216)
    ax.set_ylim(0, 145)
    ax.axis("off")
    fig.text(0.5, 0.975, "qamodel.phmm  —  objects, entry points, engines",
             ha="center", fontname=SERIF, fontsize=15, fontweight="bold",
             color=INK)
    fig.text(0.5, 0.949,
             "value objects in, value objects out; two engines that must agree",
             ha="center", fontname=SERIF, fontsize=10, color="#5d646e")

    # ---- tier 1: what you pass in -------------------------------------
    box(ax, 6, 116, 52, 24, "params.py", [
        "PhmmParams   frozen, 5 floats",
        "  delta_i delta_d gamma_i gamma_d tau",
        "  .log_trans() -> LogTrans",
        "LogTrans     11 log edges, no b_d/d_e",
        "QualityTable log_match[94] log_mismatch[94]",
        "  .cost() -> s per Phred",
        "quality_table(log_insert=None)",
    ], "in")

    box(ax, 68, 116, 54, 24, "band.py", [
        "Band         lo[L+1], hi[L+1]",
        "  .diagonal(L,T,olo,ohi,width=100)",
        "  .from_cigar(L,T,cigar,ref_start,w)",
        "  .max_width  .rows  .covers(u,v)",
        "cigar_path(L, cigar, ref_start)",
        "             -> v[u], length L+1",
    ], "in")

    box(ax, 132, 116, 74, 24, "caller supplies", [
        "read      str | uint8[]   as sequenced",
        "qual      int | int[]     Phred, same order",
        "template  str | uint8[]   from ONE origin,",
        "          strand-oriented, length ~L + slack",
        "origin    None | int | (lo, hi)",
        "          which start positions this call owns",
    ], "in", title_mono=False)

    # ---- tier 2: public API, grouped by the engine each one calls ------
    api = [
        (8, "forward()", ["oracle: full table", "O(L*T) memory",
                          "keeps full_template", "small inputs only"]),
        (72, "viterbi()", ["best single path", "with traceback",
                           "never used to score", "locates the read"]),
        (136, "forward_window()", ["two rows + banded", "O(T) memory",
                                   "no full_template: it is", "wrong for scoring"]),
    ]
    api_pts = {}
    for x, name, lines in api:
        api_pts[name] = box(ax, x, 76, 58, 22, name, lines, "api")

    # ---- tier 3: engines ----------------------------------------------
    box(ax, 8, 34, 84, 32, "_dp_full_table(...)", [
        "fM, fI, fD  shape (L+1, T+1)   + optional pointers",
        "use_max   False -> log-sum-exp   True -> max",
        "want_ptr  back-pointers for traceback",
        "",
        "the obvious statement of the recursion; the scalar D",
        "loop is kept deliberately, so it stays the oracle",
        "",
        "helpers:  _log_add()  _terminal_row()  _score_from_terminal()",
    ], "eng")

    box(ax, 112, 34, 92, 32, "_dp_two_rows(...)", [
        "prev_m/i/d  <->  cur_m/i/d, swapped each read base",
        "row u reads only row u-1, so row u-2 is dead",
        "5.1 GB -> 9 KB for a 15 kb read",
        "",
        "band=None -> every column;  band -> lo[u]..hi[u]",
        "returns (terminal_row, band_edge)",
        "",
        "helpers:  _deletion_row()  _gather()  _logsumexp()",
    ], "eng")

    # ---- tier 4: outputs ----------------------------------------------
    box(ax, 8, 6, 72, 20, "ViterbiResult", [
        "logp     float   best single path",
        "origin   int     where it started",
        "end      int     last template base used",
        "path     str     'MMMMDMMM...'",
    ], "out")

    box(ax, 100, 6, 104, 20, "WindowScore", [
        "logp        float   log P(read | this position)",
        "terminal    array   mass by END position",
        "band_edge   float   worst row's edge share, 0 if unbanded",
        ".edge_mass(margin)   read running off the right",
    ], "out")

    # ---- wiring --------------------------------------------------------
    arrow(ax, (32, 116), (37, 99))
    arrow(ax, (163, 116), (101, 99), rad=0.10, label="read, qual,\ntemplate, origin", lp=0.55)
    arrow(ax, (95, 116), (158, 99), rad=-0.16, label="band=", lp=0.66)

    arrow(ax, (37, 76), (34, 66.5), label="_prepare", lp=0.55)
    arrow(ax, (101, 76), (66, 66.5), label="use_max=True\nwant_ptr=True", lp=0.55)
    arrow(ax, (165, 76), (150, 66.5), label="_prepare", lp=0.55)

    arrow(ax, (34, 34), (40, 26.5))
    arrow(ax, (74, 34), (118, 26.5), rad=-0.06)
    arrow(ax, (150, 34), (150, 26.5))

    ax.add_patch(FancyArrowPatch((92, 50), (112, 50), arrowstyle="<|-|>",
                                 mutation_scale=12, linewidth=1.2,
                                 color="#b0563f", linestyle=(0, (5, 3)), zorder=1))
    ax.text(102, 54.5, "must agree\nto 1e-9", ha="center", va="center",
            fontname=SERIF, fontsize=8.4, color="#b0563f", zorder=4,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="none", alpha=0.95))

    ax.text(108, 1.0,
            "Everything above imports nothing from qatool: the dependency runs "
            "one way only, qatool.plotters.likelihood -> qamodel.",
            ha="center", va="center", fontname=SERIF, fontsize=9,
            color="#5d646e")

    fig.savefig(path, dpi=140)
    print(f"wrote {path}")


if __name__ == "__main__":
    draw(sys.argv[1] if len(sys.argv) > 1 else "docs/phmm_architecture.png")
