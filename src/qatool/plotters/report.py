"""Combined report figure: k-mer dotplot background, read-support overlay,
contig color bars, and marginal coverage plots.

`report_plot` orchestrates the panels; the private `_draw_*` / `_load_*` helpers
own one phase each so the figure assembly stays readable.
"""

import math
import time
import sqlite3

import pysam
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.transforms import blended_transform_factory, offset_copy
from mpl_toolkits.axes_grid1 import make_axes_locatable

from qatool.SegmentStorage import SegmentCollection
from qatool.plotters.colors import (
    SUPPORT_COLORS, OCCURRENCE_COLORS,
    PID_TIE_EPS_DEFAULT, BEST_PID_DELTA_DEFAULT,
)
from qatool.plotters.axis import _build_axis_map, _coord_ticks
from qatool.plotters.db import _fetch_aln_from_db, _collect_bam_read_names
from qatool.plotters.classify import (
    _classify_reads_for_support_colors,
    _classify_reads_for_occurrence_colors,
)
from qatool.plotters.kmer import dotplot_kmer_mosaic
from qatool.plotters.coverage import alignment_plot
from qatool.plotters.reads_overlay import reads_alignment_dotplot
from qatool.plotters.likelihood import likelihood_report


def _resolve_binsize(seq_sc, asm1: str, asm2: str, k: int, binsize: int | None) -> int:
    """Auto-compute binsize (~1000 bins on the longer axis, never below k) if unset."""
    if binsize is not None:
        return binsize
    total_a = sum(seg.end - seg.start for seg in seq_sc if seg.assembly == asm1)
    total_b = sum(seg.end - seg.start for seg in seq_sc if seg.assembly == asm2)
    return max(k, math.ceil(max(total_a, total_b, 1) / 1000))


def _load_sequences(seq_sc, asm1: str, asm2: str,
                    asm1_fasta: str | None, asm2_fasta: str | None,
                    seq_asm1: str, seq_asm2: str) -> tuple[str, str]:
    """Fetch per-assembly concatenated sequences from FASTA in seq_sc order.

    When fasta paths are given they take precedence and the passed seq_asm1/seq_asm2
    are ignored — this guarantees sequence order always matches the axis map.
    """
    if asm1_fasta is None or asm2_fasta is None:
        return seq_asm1, seq_asm2
    fa1 = pysam.FastaFile(asm1_fasta)
    fa2 = pysam.FastaFile(asm2_fasta)
    try:
        p1, p2 = [], []
        for seg in seq_sc:
            if seg.assembly == asm1:
                p1.append(fa1.fetch(seg.contig, seg.start, seg.end))
            elif seg.assembly == asm2:
                p2.append(fa2.fetch(seg.contig, seg.start, seg.end))
        return "".join(p1), "".join(p2)
    finally:
        fa1.close()
        fa2.close()


def _split_sequence(seq: str, segs_info: list) -> list:
    """Split a concatenated sequence into per-segment strings (raw lengths)."""
    parts, off = [], 0
    for _, start, end, _ in segs_info:
        parts.append(seq[off: off + (end - start)])
        off += end - start
    return parts


def _collect_all_reads(read_names, asm2_bam, asm1_bam,
                       asm2_segs_list, asm1_segs_list,
                       pid_threshold, plen_threshold) -> set:
    """Expand the seed read set with all BAM reads touching the marginal segments."""
    all_reads = set(read_names)
    if asm2_bam is not None:
        all_reads |= _collect_bam_read_names(asm2_bam, asm2_segs_list, pid_threshold, plen_threshold)
    if asm1_bam is not None:
        all_reads |= _collect_bam_read_names(asm1_bam, asm1_segs_list, pid_threshold, plen_threshold)
    return all_reads


def _classify_all(db_path, all_reads, pid_threshold, plen_threshold,
                  asm1, asm2, pid_tie_eps, aln_by_asm,
                  asm1_fasta, asm2_fasta, best_pid_delta) -> dict:
    """Support-classify every read globally (segs_by_asm=None) so subset views keep
    each read's globally-correct category.  Opens the FASTAs for the identity check."""
    fa1 = pysam.FastaFile(asm1_fasta) if asm1_fasta else None
    fa2 = pysam.FastaFile(asm2_fasta) if asm2_fasta else None
    try:
        return _classify_reads_for_support_colors(
            db_path, all_reads, pid_threshold, plen_threshold,
            segs_by_asm=None, asm1=asm1, asm2=asm2, pid_tie_eps=pid_tie_eps,
            aln_by_asm=aln_by_asm, asm1_fasta=fa1, asm2_fasta=fa2,
            best_pid_delta=best_pid_delta,
        )
    finally:
        if fa1:
            fa1.close()
        if fa2:
            fa2.close()


def _draw_contig_bars(fig, ax, divider, asm1_segs_info, asm2_segs_info,
                      plot_w, plot_h, dpi, height):
    """Colored per-contig bars below (asm2) and left (asm1) of the dotplot.
    Same-contig segments merge into one bar; bars too small to label go to a legend.
    Returns (ax_xbar, ax_ybar)."""
    all_contigs = sorted({c for c, _, _, _ in asm1_segs_info + asm2_segs_info})
    _cmap = plt.get_cmap('tab10' if len(all_contigs) <= 10 else 'tab20')
    contig_color = {c: _cmap(i % _cmap.N) for i, c in enumerate(all_contigs)}

    ax_xbar = divider.append_axes("bottom", size=0.15, pad=0.0, sharex=ax)
    ax_ybar = divider.append_axes("left",   size=0.15, pad=0.0, sharey=ax)

    # Hide everything on both bar axes except the outer edge that will carry tick marks.
    ax_xbar.set_facecolor("none")
    for sp in ("top", "left", "right"):
        ax_xbar.spines[sp].set_visible(False)
    ax_xbar.yaxis.set_visible(False)
    ax_xbar.xaxis.set_ticks_position('bottom')  # force outermost edge; sharex default is 'top'
    ax_xbar.tick_params(axis="x", bottom=True, top=False, labelbottom=False, length=4)

    ax_ybar.set_facecolor("none")
    for sp in ("top", "right", "bottom"):
        ax_ybar.spines[sp].set_visible(False)
    ax_ybar.xaxis.set_visible(False)
    ax_ybar.tick_params(axis="y", left=True, right=False, labelleft=False, length=4)

    # Caption fitting: estimate pixels per data-unit to decide if text fits in the bar.
    _bar_fontsize = 7
    _char_px = _bar_fontsize / 72 * dpi * 0.62   # ~pixels per character at this dpi
    _fig_w_px = 16.0 * dpi                        # approximate full figure width
    _fig_h_px = (height - 2.15) * dpi             # approximate main-axes height (minus bar+top)
    _legend_entries = {}                          # contig -> color, for bars too small to label

    def _contig_runs(segs_info, total):
        """Yield (contig, bar_off, bar_end) merging consecutive same-contig segments."""
        if not segs_info:
            return
        cur_contig, _, _, cur_off = segs_info[0]
        for i in range(1, len(segs_info)):
            if segs_info[i][0] != cur_contig:
                yield cur_contig, cur_off, segs_info[i][3]
                cur_contig, _, _, cur_off = segs_info[i]
        yield cur_contig, cur_off, total

    ax_xbar.set_xlim(0, plot_w)
    ax_xbar.set_ylim(0, 1)
    for contig, bar_off, bar_end in _contig_runs(asm2_segs_info, plot_w):
        bar_w = bar_end - bar_off
        ax_xbar.add_patch(plt.Rectangle(
            (bar_off, 0), bar_w, 1,
            color=contig_color[contig], transform=ax_xbar.transData, clip_on=True
        ))
        if bar_w / plot_w * _fig_w_px >= len(contig) * _char_px * 1.4:
            ax_xbar.text(bar_off + bar_w / 2, 0.5, contig,
                         ha='center', va='center', fontsize=_bar_fontsize,
                         color='black', clip_on=True)
        else:
            _legend_entries[contig] = contig_color[contig]

    ax_ybar.set_xlim(0, 1)
    ax_ybar.set_ylim(0, plot_h)
    for contig, bar_off, bar_end in _contig_runs(asm1_segs_info, plot_h):
        bar_h = bar_end - bar_off
        ax_ybar.add_patch(plt.Rectangle(
            (0, bar_off), 1, bar_h,
            color=contig_color[contig], transform=ax_ybar.transData, clip_on=True
        ))
        if bar_h / plot_h * _fig_h_px >= len(contig) * _char_px * 1.4:
            ax_ybar.text(0.5, bar_off + bar_h / 2, contig,
                         ha='center', va='center', fontsize=_bar_fontsize,
                         rotation=90, color='black', clip_on=True)
        else:
            _legend_entries[contig] = contig_color[contig]

    if _legend_entries:
        _handles = [Patch(color=c, label=n) for n, c in _legend_entries.items()]
        fig.legend(handles=_handles, loc='lower center',
                   ncol=min(len(_legend_entries), 6), fontsize=8,
                   bbox_to_anchor=(0.5, 0), frameon=True)

    return ax_xbar, ax_ybar


def _draw_marginals(fig, ax, divider, db_path, read_classes, aln_by_asm,
                    asm1_segs_list, asm2_segs_list, asm2_bam, asm1_bam,
                    asm1, asm2, pid_threshold, plen_threshold, binsize,
                    best_pid_delta, classification, plot_w, plot_h,
                    max_coverage, support_stripes, verbose):
    """Marginal coverage plots (top=asm2, right=asm1), the support/occurrence
    stripes on the main dotplot, and shared coverage ticks.  Returns (ax_top, ax_right),
    or (None, None) when either BAM is missing."""
    if asm2_bam is None or asm1_bam is None:
        return None, None

    ax_top = divider.append_axes("top", size=2, pad=0, sharex=ax)
    ax_right = divider.append_axes("right", size=2, pad=0, sharey=ax)
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_right.tick_params(axis="y", labelleft=False)

    if verbose:
        likelihood_report(db_path, read_classes, asm1_segs_list, asm2_segs_list,
                          asm1, asm2, pid_threshold, plen_threshold)

    if classification == "occurrence":
        _common = read_classes.get("common", set())
        asm2_buckets = _classify_reads_for_occurrence_colors(aln_by_asm, asm2, _common, best_pid_delta=best_pid_delta)
        asm2_colors  = OCCURRENCE_COLORS
        asm1_buckets = _classify_reads_for_occurrence_colors(aln_by_asm, asm1, _common, best_pid_delta=best_pid_delta)
        asm1_colors  = OCCURRENCE_COLORS
    else:  # "support"
        _bucket_keys   = ["asm1_only", "asm2_only", "asm1_best", "asm2_best", "common", "other"]
        _bucket_colors = [SUPPORT_COLORS[k] for k in _bucket_keys]
        asm2_buckets   = [read_classes[k] for k in _bucket_keys]
        asm2_colors    = _bucket_colors
        asm1_buckets   = [read_classes[k] for k in _bucket_keys]
        asm1_colors    = _bucket_colors

    ax_top, maxx, asm2_bucket_bins = alignment_plot(
        asm2_bam, asm2_segs_list, pid_threshold=pid_threshold,
        plen_threshold=plen_threshold, binsize=binsize,
        ax=ax_top, read_buckets=asm2_buckets, read_bucket_colors=asm2_colors,
        orientation="x", return_max=True, best_pid_delta=best_pid_delta)

    ax_right, maxy, asm1_bucket_bins = alignment_plot(
        asm1_bam, asm1_segs_list, pid_threshold=pid_threshold,
        plen_threshold=plen_threshold, binsize=binsize,
        ax=ax_right, read_buckets=asm1_buckets, read_bucket_colors=asm1_colors,
        orientation="y", return_max=True, best_pid_delta=best_pid_delta)

    # Transparent stripes on the main dotplot: red vertical = asm2_only bins,
    # blue horizontal = asm1_only bins.  bucket index 1 = asm2_only, index 0 = asm1_only.
    if support_stripes and classification == "support":
        if asm2_bucket_bins is not None:
            patches = [Rectangle((i * binsize, 0), binsize, plot_h)
                       for i, v in enumerate(asm2_bucket_bins[1]) if v > 0]
            if patches:
                ax.add_collection(PatchCollection(
                    patches, facecolor=SUPPORT_COLORS["asm2_only"], alpha=0.3, linewidth=0,
                ))
        if asm1_bucket_bins is not None:
            patches = [Rectangle((0, i * binsize), plot_w, binsize)
                       for i, v in enumerate(asm1_bucket_bins[0]) if v > 0]
            if patches:
                ax.add_collection(PatchCollection(
                    patches, facecolor=SUPPORT_COLORS["asm1_only"], alpha=0.3, linewidth=0,
                ))

    m = max_coverage if max_coverage is not None else max(maxx, maxy)
    ax_top.set_ylim(0, m)
    ax_right.set_xlim(0, m)
    # Shared tick positions so both coverage axes are identical.
    cov_ticks = mticker.MaxNLocator(nbins=5, integer=True).tick_values(0, m)
    cov_ticks = [t for t in cov_ticks if 0 <= t <= m]
    ax_top.set_yticks(cov_ticks)
    ax_right.set_xticks(cov_ticks)
    ax_top.set_xticks([])
    ax_top.set_xlabel("")
    ax_right.set_yticks([])
    ax_right.set_ylabel("")

    return ax_top, ax_right


def _draw_coordinate_axes(fig, ax, ax_xbar, ax_ybar, ax_top, ax_right,
                          asm1_segs_info, asm2_segs_info, plot_w, plot_h):
    """Genome-coordinate ticks and per-segment labels on the color-bar axes.

    Set LAST: ax_xbar/ax_ybar share the major-tick container with ax/ax_top (sharex)
    and ax/ax_right (sharey), so any earlier set_xticks([]) would wipe these.
    """
    ax.set_xlabel("")
    ax.set_ylabel("")

    x_major, x_minor = _coord_ticks(asm2_segs_info, plot_w)
    ax_xbar.set_xticks(x_major)
    ax_xbar.xaxis.set_minor_locator(mticker.FixedLocator(x_minor))
    ax_xbar.tick_params(axis="x", which="major", bottom=True, top=False, labelbottom=False, length=5)
    ax_xbar.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, length=3)
    # Segment labels at the center of each segment, shifted down a bit below the bar.
    _xtrans = offset_copy(blended_transform_factory(ax_xbar.transData, ax_xbar.transAxes),
                          fig=fig, x=0, y=-4, units='points')
    for _, gen_s, gen_e, off in asm2_segs_info:
        ax_xbar.text(off + (gen_e - gen_s) / 2, 0, f"{gen_s:,}:{gen_e:,}",
                     transform=_xtrans, ha='center', va='top',
                     fontsize=7, rotation=45, clip_on=False)
    # Suppress both major and minor on all other x-sharing axes.
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=False)
    if ax_top is not None:
        ax_top.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=False)

    y_major, y_minor = _coord_ticks(asm1_segs_info, plot_h)
    ax_ybar.set_yticks(y_major)
    ax_ybar.yaxis.set_minor_locator(mticker.FixedLocator(y_minor))
    ax_ybar.tick_params(axis="y", which="major", left=True, right=False, labelleft=False, length=5)
    ax_ybar.tick_params(axis="y", which="minor", left=True, right=False, labelleft=False, length=3)
    # Segment labels at the center of each segment, shifted right toward the bar.
    _ytrans = offset_copy(blended_transform_factory(ax_ybar.transAxes, ax_ybar.transData),
                          fig=fig, x=-8, y=0, units='points')
    for _, gen_s, gen_e, off in asm1_segs_info:
        ax_ybar.text(0, off + (gen_e - gen_s) / 2, f"{gen_s:,}:{gen_e:,}",
                     transform=_ytrans, ha='right', va='center',
                     fontsize=7, clip_on=False)
    # Suppress both major and minor on all other y-sharing axes.
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
    if ax_right is not None:
        ax_right.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)


def report_plot(seq_asm1: str, seq_asm2: str, segments: SegmentCollection, read_names: set[str],
                db_path: str, pid_threshold: float, plen_threshold: float, k: int = 15,
                binsize: int | None = None,
                output_path: str = "report_plot.png", dpi: int = 150,
                asm2_bam: str | None = None, asm1_bam: str | None = None,
                asm1_fasta: str | None = None, asm2_fasta: str | None = None,
                boundary_params: tuple[str, float, float] = ("black", 1.0, 1.0),
                asm1: str = "asm1", asm2: str = "asm2",
                pid_tie_eps: float = PID_TIE_EPS_DEFAULT,
                best_pid_delta: float = BEST_PID_DELTA_DEFAULT,
                max_kmer_freq: int = 100,
                max_coverage: float | None = None, show: bool = True,
                support_stripes: bool = True,
                seq_segments: SegmentCollection | None = None,
                title: str | None = None,
                classification: str = "support",
                verbose: bool = False):
    """
    Create a combined report figure (k-mer dotplot background + read-alignment
    overlay + contig bars + marginal coverage) and save it to output_path.

    Parameters (selected):
    - max_coverage: fix the coverage axis limit for both marginal alignment plots
      (ax_top y-axis and ax_right x-axis).  When None (default), the limit is set to
      max(maxx, maxy) from the data — shared within a call but possibly differing
      across calls.  Pass the same value to all calls to make plots comparable.
    - show: when True, display the figure with plt.show(); set False for batch runs.
    - verbose: when True, print per-step timings and the likelihood report; the report
      is silent by default.  Run likelihood_report() directly to get the numbers.

    Returns the matplotlib Figure.
    """
    def _t(label, t0):
        if verbose:
            print(f"  [report_plot] {label}: {time.perf_counter()-t0:.2f}s")
        return time.perf_counter()

    _t0 = time.perf_counter()

    sc_sorted = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    # seq_sc defines the coordinate system for the sequences, the dotplot extent,
    # alignment-line mapping, and contig bars.  When seq_segments is provided (sequences
    # were built from different segments than `segments`), use it so the alignment lines
    # stay consistent with the dotplot content; `segments` is kept for marginal BAM queries.
    seq_sc = (sorted(seq_segments, key=lambda s: (s.contig, s.start, s.end))
              if seq_segments is not None else sc_sorted)

    binsize = _resolve_binsize(seq_sc, asm1, asm2, k, binsize)

    # Bin-aligned axis maps: each segment padded to a whole number of bins so segment
    # boundaries are exact multiples of binsize.
    asm1_segs_info, plot_h = _build_axis_map(seq_sc, asm1, binsize=binsize)
    asm2_segs_info, plot_w = _build_axis_map(seq_sc, asm2, binsize=binsize)

    # Segment lists used by marginal plots and BAM read collection (always from segments).
    asm2_segs_list = [s for s in sc_sorted if s.assembly == asm2]
    asm1_segs_list = [s for s in sc_sorted if s.assembly == asm1]

    seq_asm1, seq_asm2 = _load_sequences(seq_sc, asm1, asm2, asm1_fasta, asm2_fasta,
                                         seq_asm1, seq_asm2)

    _t0 = _t(f"setup  seqs={len(seq_asm1)//1000}kbp/{len(seq_asm2)//1000}kbp  "
             f"reads={len(read_names)}  segs={len(sc_sorted)}", _t0)

    # Expand read set so classification covers every read visible in the marginals.
    all_reads = _collect_all_reads(read_names, asm2_bam, asm1_bam,
                                   asm2_segs_list, asm1_segs_list,
                                   pid_threshold, plen_threshold)
    _t0 = _t(f"collect_bam_read_names  ({len(all_reads) - len(read_names):+d} reads)", _t0)

    # Fetch alignments once for all reads; reused by dotplot overlay and classifier.
    con = sqlite3.connect(db_path)
    try:
        aln_by_asm = _fetch_aln_from_db(con, list(all_reads), [asm1, asm2], pid_threshold, plen_threshold)
    finally:
        con.close()
    _t0 = _t(f"fetch_aln_from_db  ({sum(len(v) for v in aln_by_asm.values())} rows)", _t0)

    # Classify reads early so the overlay and marginals share the same read_classes.
    read_classes = None
    if asm2_bam is not None and asm1_bam is not None:
        read_classes = _classify_all(db_path, all_reads, pid_threshold, plen_threshold,
                                     asm1, asm2, pid_tie_eps, aln_by_asm,
                                     asm1_fasta, asm2_fasta, best_pid_delta)
        _t0 = _t("classify_reads_for_support_colors", _t0)

    # Split concatenated sequences into per-segment lists for the mosaic dotplot.
    seqs_a = _split_sequence(seq_asm1, asm1_segs_info)
    seqs_b = _split_sequence(seq_asm2, asm2_segs_info)

    # Figure height: main axes is bins_a/bins_b tall for square dotplot pixels, +2in for ax_top.
    total_bins_a = plot_h // binsize
    total_bins_b = plot_w // binsize
    height = max(4.0, 16.0 * total_bins_a / max(total_bins_b, 1) + 2)
    fig, ax = plt.subplots(figsize=(16.0, height), dpi=dpi)

    # Background k-mer mosaic dotplot (no mixed bins at segment boundaries).
    dotplot_kmer_mosaic(seqs_a, seqs_b, k=k, binsize=binsize, ax=ax,
                        max_kmer_freq=max_kmer_freq, verbose=verbose)
    ax.set_aspect('auto')
    _t0 = _t("dotplot_kmer_mosaic", _t0)

    # Full-size transparent inset for the reads overlay.
    ax_overlay = ax.inset_axes((0, 0, 1, 1))
    ax_overlay.set_facecolor("none")
    ax_overlay.patch.set_alpha(0.0)
    ax_overlay.set_frame_on(False)
    b_color, b_alpha, b_lw = boundary_params
    reads_alignment_dotplot(
        seq_sc, all_reads, pid_threshold, plen_threshold, db_path,
        ax=ax_overlay, boundary_color=b_color, boundary_alpha=b_alpha, boundary_lw=b_lw,
        asm1=asm1, asm2=asm2, pid_tie_eps=pid_tie_eps, aln_by_asm=aln_by_asm,
        binsize=binsize, read_classes=read_classes,
    )
    _t0 = _t("reads_alignment_dotplot", _t0)

    # Keep overlays aligned with base axes.
    ax.set_xlim(0, plot_w)
    ax.set_ylim(0, plot_h)
    ax_overlay.set_xlim(ax.get_xlim())
    ax_overlay.set_ylim(ax.get_ylim())

    divider = make_axes_locatable(ax)
    ax_xbar, ax_ybar = _draw_contig_bars(fig, ax, divider, asm1_segs_info, asm2_segs_info,
                                         plot_w, plot_h, dpi, height)

    ax_top, ax_right = _draw_marginals(fig, ax, divider, db_path, read_classes, aln_by_asm,
                                       asm1_segs_list, asm2_segs_list, asm2_bam, asm1_bam,
                                       asm1, asm2, pid_threshold, plen_threshold, binsize,
                                       best_pid_delta, classification, plot_w, plot_h,
                                       max_coverage, support_stripes, verbose)
    if ax_top is not None:
        _t0 = _t("draw_marginals", _t0)

    # Clear ticks and labels from the main dotplot's transparent overlay.
    ax_overlay.set_xticks([])
    ax_overlay.set_yticks([])
    ax_overlay.xaxis.set_minor_locator(mticker.NullLocator())
    ax_overlay.yaxis.set_minor_locator(mticker.NullLocator())
    ax_overlay.set_xlabel("")
    ax_overlay.set_ylabel("")

    _draw_coordinate_axes(fig, ax, ax_xbar, ax_ybar, ax_top, ax_right,
                          asm1_segs_info, asm2_segs_info, plot_w, plot_h)

    if title:
        fig.suptitle(title, fontsize=10, y=1.01)
    fig.savefig(output_path, bbox_inches="tight")
    _t0 = _t("savefig", _t0)
    if show:
        plt.show()
    plt.close(fig)
    return fig
