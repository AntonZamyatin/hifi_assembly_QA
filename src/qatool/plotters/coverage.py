"""Marginal coverage / alignment bar plot along concatenated segments."""

import math

import numpy as np
import pysam
import matplotlib.pyplot as plt

from qatool.aln_utils import cigar_stats
from qatool.SegmentStorage import Segment
from qatool.plotters.colors import BEST_PID_DELTA_DEFAULT


def alignment_plot(bamfile: str, segments: list,
                   pid_threshold: float = 0.99, plen_threshold: float | None = None,
                   binsize: int = 1000,
                   ax=None, read_buckets=None, read_bucket_colors=None, orientation: str = "x",
                   return_max=False, best_pid_delta: float = BEST_PID_DELTA_DEFAULT):
    """
    Generate a simple alignment plot for a given BAM file, contig, and region.

    Parameters:
    - bamfile: Path to the BAM file containing alignments.
    - segments: List of Segment to plot alignments for (same assembly/contig).
    """
    if not segments:
        raise ValueError("segments must be a non-empty list")
    segments = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    # Each segment occupies ceil(seg_len/binsize) bins so that segment offsets
    # land exactly on bin boundaries (no bin straddles a segment edge).
    bins = [0.0 for _ in range(sum(math.ceil((s.end - s.start) / binsize) for s in segments))]
    bucket_bins = None
    name_to_bucket = None
    if read_buckets is not None:
        bucket_bins = [[0.0 for _ in range(len(bins))] for _ in range(len(read_buckets))]
        name_to_bucket = {}
        for i, bucket in enumerate(read_buckets):
            for name in bucket:
                name_to_bucket[name] = i

    # Collect valid alignments; later keep only best-PID per read.
    valid_alns = []
    offset = 0
    bam = pysam.AlignmentFile(bamfile, "rb")
    try:
        for seg in segments:
            for aln in bam.fetch(seg.contig, start=seg.start, end=seg.end):
                soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid, _ = cigar_stats(aln)
                if pid < pid_threshold:
                    continue
                if plen_threshold is not None and (plen is None or plen < plen_threshold):
                    continue
                startbin = math.floor((max(aln.reference_start, seg.start) - seg.start + offset) / binsize)
                endbin = math.ceil((min(aln.reference_end, seg.end) - seg.start + offset) / binsize)
                if startbin < endbin:
                    valid_alns.append((aln.query_name, float(pid), startbin, endbin))
            offset += math.ceil((seg.end - seg.start) / binsize) * binsize
    finally:
        bam.close()

    # Pass 1: find best pid per read.
    best_pid_by_read = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if read_name not in best_pid_by_read or pid > best_pid_by_read[read_name]:
            best_pid_by_read[read_name] = pid

    # Pass 2: collect qualifying alignments (within best_pid_delta of best) per read,
    # then fill bins with weight 1 / number of tied alignments.
    qualifying: dict = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if pid + best_pid_delta >= best_pid_by_read.get(read_name, float("-inf")):
            qualifying.setdefault(read_name, []).append((startbin, endbin))

    for read_name, ranges in qualifying.items():
        w = 1.0 / len(ranges)
        if bucket_bins is None:
            for startbin, endbin in ranges:
                for i in range(startbin, endbin):
                    bins[i] += w
        else:
            bucket_idx = name_to_bucket.get(read_name)
            target = bucket_bins[bucket_idx] if bucket_idx is not None else bins
            for startbin, endbin in ranges:
                for i in range(startbin, endbin):
                    target[i] += w

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(16, 5))

    # n+1 x-values so step='post' renders all n bins fully:
    # the last bin covers [x[n-1], x[n]] and needs a right edge.
    x = np.arange(len(bins) + 1) * binsize
    if bucket_bins is None:
        bins_plot = bins + [0.0]
        if orientation == "x":
            ax.stackplot(x, bins_plot, color='blue', step='post')
        elif orientation == "y":
            ax.fill_betweenx(x, 0, bins_plot, step='post', color='blue', linewidth=0)
        else:
            raise ValueError("orientation must be 'x' or 'y'")
    else:
        if read_bucket_colors is None:
            colors = [f"C{i % 10}" for i in range(len(bucket_bins))] + ['#9e9e9e']
        else:
            if len(read_bucket_colors) != len(bucket_bins):
                raise ValueError("read_bucket_colors length must match read_buckets length")
            colors = list(read_bucket_colors) + ['#6e6e6e']
        series = [s + [0.0] for s in bucket_bins + [bins]]
        if orientation == "x":
            ax.stackplot(x, series, colors=colors, step='post')
        elif orientation == "y":
            y0 = np.zeros(len(x))
            for s, c in zip(series, colors):
                y1 = y0 + np.array(s)
                ax.fill_betweenx(x, y0, y1, step='post', color=c, linewidth=0)
                y0 = y1
        else:
            raise ValueError("orientation must be 'x' or 'y'")

    # Axes styling: no rectangle, axes through origin
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    if standalone:
        if orientation == "x":
            ax.set_xlim(0, max(len(bins) * binsize - 1, 0))
        else:
            ax.set_ylim(0, max(len(bins) * binsize - 1, 0))
        plt.xlabel('Position in Bin')
        plt.ylabel('Number of Alignments')
        plt.title('Alignment Plot for Segments (' + str(len(segments)) + ')')
        plt.show()

    if return_max:
        if bucket_bins is None:
            return ax, max(bins) if bins else 0.0, None
        total = np.array(bins, dtype=float)
        for s in bucket_bins:
            total += np.array(s, dtype=float)
        return ax, float(np.max(total)) if len(total) > 0 else 0.0, bucket_bins
    return ax
