"""
Project: hifi assembly QA
Module: alignment_plot.py
Description: Dot-plot style visualization of read alignments between segment collections.

Author: Anton Zamyatin
Affiliation: Pennsylvania State University
Contact: zamyatin.anton@gmail.com

Created: 2/12/26

License: MIT License
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
import matplotlib.pyplot as plt
import sqlite3
import pysam
import math

from qatool.SegmentStorage import Segment, SegmentCollection
from qatool.aln_utils import *


def simple_alignment_plot(bamfile: str, contig: str, start: int | None = None,
                          end: int | None = None, pid_threshold: float = 0.99,
                          binsize: int = 1000):
    """
    Generate a simple alignment plot for a given BAM file, contig, and region.

    Parameters:
    - bamfile: Path to the BAM file containing alignments.
    - contig: Contig name to plot alignments for.
    - start: Start position of the region to plot.
    - end: End position of the region to plot.
    """
    bam = pysam.AlignmentFile(bamfile, "rb")
    alns = bam.fetch(contig, start=start, end=end)
    if not start: start = 0
    if not end: end = bam.get_reference_length(contig)
    bins = [0 for i in range(math.ceil((end-start)/binsize))]
    for aln in alns:
        pid = cigar_stats(aln)[-2]
        if pid >= pid_threshold:
            startbin = math.floor((aln.reference_start-start)/binsize)
            endbin = math.ceil((aln.reference_end-start)/binsize)
            for i in range(startbin, endbin): bins[i] += 1
    fig = plt.figure(figsize=(16, 5))
    plt.stackplot(range(len(bins)), bins, color='blue', step='pre')
    plt.xlabel('Position in Bin')
    plt.ylabel('Number of Alignments')
    plt.title('Alignment Plot for ' + contig + ' from ' + str(start) + ' to ' + str(end))
    plt.show()





@dataclass(frozen=True)
class AlnRecord:
    assembly: str
    read_name: str
    contig: str
    start: int
    end: int
    strand: Optional[int]
    aln_type: str
    pid: Optional[float]


@dataclass(frozen=True)
class AxisSegment:
    segment: Segment
    offset: int
    length: int

    def coord(self, pos: int) -> int:
        return self.offset + (pos - self.segment.start)


class SegmentAxis:
    def __init__(self, segments: Sequence[Segment], gap: int):
        ordered = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
        self.gap = gap
        self.items: List[AxisSegment] = []
        offset = 0
        for seg in ordered:
            length = len(seg)
            self.items.append(AxisSegment(seg, offset, length))
            offset += length + gap
        self.total_len = max(offset - gap, 0)
        self._by_contig: Dict[str, List[AxisSegment]] = {}
        for item in self.items:
            self._by_contig.setdefault(item.segment.contig, []).append(item)

    def find_axis_segment(self, contig: str, start: int, end: int) -> Optional[AxisSegment]:
        """Return AxisSegment if the interval is fully contained in a segment on this axis."""
        for item in self._by_contig.get(contig, []):
            seg = item.segment
            if start >= seg.start and end <= seg.end:
                return item
        return None

    def centers_and_labels(self, max_labels: int) -> Tuple[List[float], List[str]]:
        if len(self.items) > max_labels:
            return [], []
        centers = []
        labels = []
        for item in self.items:
            centers.append(item.offset + item.length / 2)
            seg = item.segment
            labels.append(f"{seg.contig}:{seg.start}-{seg.end}")
        return centers, labels

    def boundaries(self) -> List[int]:
        bounds = []
        for i, item in enumerate(self.items):
            if i == len(self.items) - 1:
                continue
            bounds.append(item.offset + item.length)
        return bounds


def _chunks(seq: Sequence[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(seq), size):
        yield list(seq[i:i + size])


def _fetch_alignments(
    con: sqlite3.Connection,
    reads: Sequence[str],
    assemblies: Sequence[str],
) -> Dict[str, List[AlnRecord]]:
    if not reads:
        return {}

    asm_placeholders = ",".join(["?"] * len(assemblies))
    base_sql = f"""
        SELECT assembly, read_name, contig, start, end, strand, aln_type, pid
        FROM aln
        WHERE assembly IN ({asm_placeholders})
          AND contig IS NOT NULL
          AND start IS NOT NULL
          AND end IS NOT NULL
          AND read_name IN ({{READS_IN}})
    """

    out: Dict[str, List[AlnRecord]] = {asm: [] for asm in assemblies}
    for chunk in _chunks(list(reads), 800):
        reads_in = ",".join(["?"] * len(chunk))
        sql = base_sql.replace("{READS_IN}", reads_in)
        params = list(assemblies) + chunk
        rows = con.execute(sql, params).fetchall()
        for asm, read_name, contig, start, end, strand, aln_type, pid in rows:
            out[asm].append(
                AlnRecord(str(asm), str(read_name), str(contig), int(start), int(end),
                          None if strand is None else int(strand), str(aln_type),
                          None if pid is None else float(pid))
            )
    return out


def plot_alignment_dotplot(
    segments: SegmentCollection,
    reads: Set[str],
    db_path: str,
    assembly_x: str = "cand",
    assembly_y: str = "ref",
    gap: int = 0,
    max_reads: int = 2000,
    pid_thresholds: Tuple[float, float] = (0.99, 0.99), # (cand, ref)
    show_strand: bool = True,
    max_labels: int = 20,
    line_alpha: float = 0.35,
    line_width: float = 0.6,
    figsize: Tuple[float, float] = (12, 10),
    show_marginals: bool = True,
    coverage_bins: int = 200,
    coverage_alpha: float = 0.6,
    layout_tight: bool = True,
    title_pad: float = 12.0,
    use_suptitle: bool = True,
    suptitle_y: float = 0.98,
):
    """
    Dot-plot style visualization of alignments between two assemblies.

    Axes are constructed by concatenating segments from the SegmentCollection for each assembly.
    Each read contributes short line segments between its alignments on the two axes.
    Line color is based on alignment pid difference: pid_x - pid_y.

    Returns (fig, ax, stats) where stats includes counts for plotted reads and alignments.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError("matplotlib is required for plotting") from e

    segs_x = [s for s in segments if s.assembly == assembly_x]
    segs_y = [s for s in segments if s.assembly == assembly_y]
    if not segs_x or not segs_y:
        raise ValueError("segments must include at least one segment for each assembly")

    axis_x = SegmentAxis(segs_x, gap=gap)
    axis_y = SegmentAxis(segs_y, gap=gap)

    read_list = sorted(reads)
    if max_reads is not None and len(read_list) > max_reads:
        read_list = read_list[:max_reads]

    con = sqlite3.connect(db_path)
    try:
        aln_by_asm = _fetch_alignments(con, read_list, [assembly_x, assembly_y])
    finally:
        con.close()

    # Index alignments by read for fast pairing
    by_read_x: Dict[str, List[AlnRecord]] = {}
    by_read_y: Dict[str, List[AlnRecord]] = {}
    for aln in aln_by_asm.get(assembly_x, []):
        by_read_x.setdefault(aln.read_name, []).append(aln)
    for aln in aln_by_asm.get(assembly_y, []):
        by_read_y.setdefault(aln.read_name, []).append(aln)

    if show_marginals:
        fig, ax = plt.subplots(figsize=figsize, dpi=120)
        try:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
        except ImportError as e:
            raise ImportError("matplotlib is required for plotting") from e
        divider = make_axes_locatable(ax)
        ax_top = divider.append_axes("top", size="25%", pad=0.1, sharex=ax)
        ax_right = divider.append_axes("right", size="25%", pad=0.1, sharey=ax)
        ax_top.tick_params(axis="x", labelbottom=False)
        ax_right.tick_params(axis="y", labelleft=False)
    else:
        fig, ax = plt.subplots(figsize=figsize, dpi=120)
        ax_top = None
        ax_right = None

    plotted_reads = 0
    plotted_lines = 0
    line_segments: List[Tuple[List[int], List[int], float]] = []
    cov_intervals_x: List[Tuple[int, int]] = []
    cov_intervals_y: List[Tuple[int, int]] = []

    pid_thr_x, pid_thr_y = pid_thresholds
    for read_name in read_list:
        alns_x = by_read_x.get(read_name)
        alns_y = by_read_y.get(read_name)
        if not alns_x or not alns_y:
            continue

        # Filter alignments to those fully contained in the axis segments
        valid_x = []
        for aln in alns_x:
            if aln.pid is None or aln.pid < pid_thr_x:
                continue
            axis_seg = axis_x.find_axis_segment(aln.contig, aln.start, aln.end)
            if axis_seg:
                valid_x.append((aln, axis_seg))
        valid_y = []
        for aln in alns_y:
            if aln.pid is None or aln.pid < pid_thr_y:
                continue
            axis_seg = axis_y.find_axis_segment(aln.contig, aln.start, aln.end)
            if axis_seg:
                valid_y.append((aln, axis_seg))

        best_x = None
        if valid_x:
            best_x = max(valid_x, key=lambda p: p[0].pid if p[0].pid is not None else -1.0)
            bx1 = best_x[1].coord(best_x[0].start)
            bx2 = best_x[1].coord(best_x[0].end)
            if show_strand and best_x[0].strand == 0:
                bx1, bx2 = bx2, bx1
            cov_intervals_x.append((bx1, bx2))

        best_y = None
        if valid_y:
            best_y = max(valid_y, key=lambda p: p[0].pid if p[0].pid is not None else -1.0)
            by1 = best_y[1].coord(best_y[0].start)
            by2 = best_y[1].coord(best_y[0].end)
            if show_strand and best_y[0].strand == 0:
                by1, by2 = by2, by1
            cov_intervals_y.append((by1, by2))

        if not valid_x or not valid_y:
            continue

        plotted_reads += 1
        for aln_x, seg_x in valid_x:
            x1 = seg_x.coord(aln_x.start)
            x2 = seg_x.coord(aln_x.end)
            if show_strand and aln_x.strand == 0:
                x1, x2 = x2, x1
            for aln_y, seg_y in valid_y:
                y1 = seg_y.coord(aln_y.start)
                y2 = seg_y.coord(aln_y.end)
                if show_strand and aln_y.strand == 0:
                    y1, y2 = y2, y1
                pid_x = aln_x.pid if aln_x.pid is not None else 0.0
                pid_y = aln_y.pid if aln_y.pid is not None else 0.0
                line_segments.append(([x1, x2], [y1, y2], pid_x - pid_y))
                plotted_lines += 1

    if line_segments:
        try:
            from matplotlib import colors as mcolors
        except ImportError as e:
            raise ImportError("matplotlib is required for plotting") from e
        max_abs = max(abs(delta) for _, _, delta in line_segments)
        max_abs = max(max_abs, 0.01)
        norm = mcolors.Normalize(vmin=-max_abs, vmax=max_abs)
        cmap = plt.get_cmap("coolwarm")
        for xs, ys, delta in line_segments:
            ax.plot(xs, ys, color=cmap(norm(delta)), alpha=line_alpha, linewidth=line_width)

    if show_marginals and ax_top is not None and ax_right is not None:
        try:
            import numpy as np
        except ImportError as e:
            raise ImportError("numpy is required for coverage plots") from e

        def _compute_coverage(axis_len: int, intervals: List[Tuple[int, int]], bins: int):
            if axis_len <= 0 or not intervals:
                return np.array([]), np.array([])
            bins = max(int(bins), 1)
            edges = np.linspace(0, axis_len, bins + 1)
            cov = np.zeros(bins, dtype=float)
            bin_w = edges[1] - edges[0]
            for start, end in intervals:
                if end < start:
                    start, end = end, start
                start = max(0.0, min(float(axis_len), float(start)))
                end = max(0.0, min(float(axis_len), float(end)))
                if end <= start:
                    continue
                i0 = int(start / axis_len * bins)
                i1 = int(end / axis_len * bins)
                if i0 >= bins:
                    continue
                if i1 < i0:
                    i1 = i0
                if i1 == i0:
                    cov[i0] += (end - start) / bin_w
                else:
                    cov[i0] += (edges[i0 + 1] - start) / bin_w
                    for i in range(i0 + 1, min(i1, bins)):
                        cov[i] += 1.0
                    if i1 < bins:
                        cov[i1] += (end - edges[i1]) / bin_w
            centers = (edges[:-1] + edges[1:]) / 2.0
            return centers, cov

        cx, covx = _compute_coverage(axis_x.total_len, cov_intervals_x, coverage_bins)
        cy, covy = _compute_coverage(axis_y.total_len, cov_intervals_y, coverage_bins)
        ax_top.set_xlim(0, max(axis_x.total_len, 1))
        ax_right.set_ylim(0, max(axis_y.total_len, 1))
        if cx.size:
            ax_top.fill_between(cx, 0, covx, color="#444444", alpha=coverage_alpha, linewidth=0)
        if cy.size:
            ax_right.fill_betweenx(cy, 0, covy, color="#444444", alpha=coverage_alpha, linewidth=0)
        ax_top.set_ylabel("coverage")
        ax_right.set_xlabel("coverage")
        ax_top.spines["right"].set_visible(False)
        ax_top.spines["top"].set_visible(False)
        ax_right.spines["top"].set_visible(False)
        ax_right.spines["right"].set_visible(False)

    # Draw contig boundaries (single line between segments)
    for x in axis_x.boundaries():
        ax.axvline(x, color="#dddddd", linewidth=0.6, zorder=0)
    for y in axis_y.boundaries():
        ax.axhline(y, color="#dddddd", linewidth=0.6, zorder=0)

    # Tick labels at segment centers (only if not too many)
    xticks, xlabels = axis_x.centers_and_labels(max_labels)
    yticks, ylabels = axis_y.centers_and_labels(max_labels)
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=8)
    if yticks:
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=8)

    ax.set_xlim(0, max(axis_x.total_len, 1))
    ax.set_ylim(0, max(axis_y.total_len, 1))
    ax.set_aspect("equal", adjustable="box")
    if show_marginals:
        try:
            data_ratio = max(axis_y.total_len, 1) / max(axis_x.total_len, 1)
            ax.set_box_aspect(data_ratio)
        except AttributeError:
            pass
    ax.set_xlabel(f"{assembly_x} segments")
    ax.set_ylabel(f"{assembly_y} segments")
    if use_suptitle:
        fig.suptitle("Read alignment dot-plot", y=suptitle_y)
    else:
        ax.set_title("Read alignment dot-plot", pad=title_pad)

    if layout_tight:
        top = 0.92 if use_suptitle else 0.95
        fig.tight_layout(rect=(0.05, 0.05, 0.98, top))

    stats = {
        "reads_requested": len(reads),
        "reads_plotted": plotted_reads,
        "lines_plotted": plotted_lines,
        "segments_x": len(segs_x),
        "segments_y": len(segs_y),
    }

    fig.savefig("plot.png")
    return fig, ax, stats
