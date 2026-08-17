"""Marginal coverage / alignment bar plot along concatenated segments.

`compute_coverage_bins` owns the BAM iteration and binning; `alignment_plot` is a
matplotlib renderer on top of it.  Keeping them separate lets other front ends
(see `qaview`) reuse the exact same numbers.
"""

import math
from dataclasses import dataclass

import numpy as np
import pysam
import matplotlib.pyplot as plt

from qatool.aln_utils import cigar_stats
from qatool.SegmentStorage import Segment
from qatool.plotters.colors import BEST_PID_DELTA_DEFAULT


@dataclass
class CoverageBins:
    """Binned coverage along concatenated segments, plus which reads land in each bin.

    Bins are laid out per segment: segment i owns ceil(seg_len / binsize) consecutive
    bins, so a bin never straddles a segment boundary.  Segments are ordered by
    (contig, start, end) — pass the same sorted order to any axis map used alongside.

    - bins        : coverage from reads in no bucket (all reads when read_buckets is None)
    - bucket_bins : per-bucket coverage, or None when read_buckets was not given
    - bin_contig / bin_start / bin_end : genomic identity of each bin (end clipped to
      the segment end, so the last bin of a segment may be shorter than binsize)
    - bin_offsets / read_idx : CSR index of bin -> reads.  The reads in bin i are
      read_idx[bin_offsets[i]:bin_offsets[i + 1]], as indices into read_names.
    """
    bins: list
    bucket_bins: list | None
    bin_contig: list
    bin_start: list
    bin_end: list
    bin_offsets: np.ndarray
    read_idx: np.ndarray
    read_names: list

    def total(self) -> np.ndarray:
        """Summed coverage across every bucket plus the unbucketed remainder."""
        out = np.array(self.bins, dtype=float)
        for series in (self.bucket_bins or []):
            out += np.array(series, dtype=float)
        return out

    def reads_in_bin(self, i: int) -> list:
        """Read names landing in bin i."""
        return [self.read_names[j] for j in self.read_idx[self.bin_offsets[i]:self.bin_offsets[i + 1]]]

    def by_read(self) -> tuple:
        """Transpose of the bin index: local read index -> bin indices, CSR-style.

        Returns (offsets, bins) with offsets of length len(read_names) + 1.  Lets a
        viewer answer "which bins does this read occupy?" without rescanning.
        """
        n = len(self.read_names)
        counts = np.zeros(n, dtype=np.int64)
        for r in self.read_idx:
            counts[r] += 1
        offsets = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(counts, out=offsets[1:])
        cursor = offsets[:-1].copy()
        bins = np.zeros(len(self.read_idx), dtype=np.int32)
        for b in range(len(self.bins)):
            for k in range(self.bin_offsets[b], self.bin_offsets[b + 1]):
                r = self.read_idx[k]
                bins[cursor[r]] = b
                cursor[r] += 1
        return offsets, bins


def compute_coverage_bins(bamfile: str, segments: list,
                          pid_threshold: float = 0.99, plen_threshold: float | None = None,
                          binsize: int = 1000, read_buckets=None,
                          best_pid_delta: float = BEST_PID_DELTA_DEFAULT) -> CoverageBins:
    """
    Bin read coverage over `segments`, tracking which reads contribute to each bin.

    A read's qualifying alignments are those within best_pid_delta of its own best
    pid; each contributes weight 1/n so a multi-mapping read still counts once in
    total.  Reads not in any bucket accumulate into `bins`.
    """
    if not segments:
        raise ValueError("segments must be a non-empty list")
    segments = sorted(segments, key=lambda s: (s.contig, s.start, s.end))

    # Per-bin genomic identity, built in the same order the bins are indexed.
    bin_contig: list = []
    bin_start: list = []
    bin_end: list = []
    for seg in segments:
        for j in range(math.ceil((seg.end - seg.start) / binsize)):
            gs = seg.start + j * binsize
            bin_contig.append(seg.contig)
            bin_start.append(gs)
            bin_end.append(min(seg.end, gs + binsize))
    n_bins = len(bin_contig)

    bins = [0.0 for _ in range(n_bins)]
    bucket_bins = None
    name_to_bucket = None
    if read_buckets is not None:
        bucket_bins = [[0.0 for _ in range(n_bins)] for _ in range(len(read_buckets))]
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

    read_names = list(qualifying)
    read_pos = {name: i for i, name in enumerate(read_names)}
    bin_reads: list = [set() for _ in range(n_bins)]

    for read_name, ranges in qualifying.items():
        w = 1.0 / len(ranges)
        if bucket_bins is None:
            target = bins
        else:
            bucket_idx = name_to_bucket.get(read_name)
            target = bucket_bins[bucket_idx] if bucket_idx is not None else bins
        ri = read_pos[read_name]
        for startbin, endbin in ranges:
            for i in range(startbin, endbin):
                target[i] += w
                bin_reads[i].add(ri)

    bin_offsets = np.zeros(n_bins + 1, dtype=np.int64)
    np.cumsum([len(s) for s in bin_reads], out=bin_offsets[1:])
    read_idx = np.fromiter((ri for s in bin_reads for ri in sorted(s)),
                           dtype=np.int32, count=int(bin_offsets[-1]))

    return CoverageBins(bins=bins, bucket_bins=bucket_bins,
                        bin_contig=bin_contig, bin_start=bin_start, bin_end=bin_end,
                        bin_offsets=bin_offsets, read_idx=read_idx, read_names=read_names)


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
    cov = compute_coverage_bins(bamfile, segments, pid_threshold=pid_threshold,
                                plen_threshold=plen_threshold, binsize=binsize,
                                read_buckets=read_buckets, best_pid_delta=best_pid_delta)
    bins = cov.bins
    bucket_bins = cov.bucket_bins

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
        total = cov.total()
        return ax, float(np.max(total)) if len(total) > 0 else 0.0, bucket_bins
    return ax
