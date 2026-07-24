"""Coordinate/axis helpers: map concatenated segments to a bin-aligned plot space."""

import math


def _build_axis_map(sc, assembly: str, binsize: int = None):
    """
    Build a (contig, start, end, offset) list for segments of the given assembly,
    in the order they appear in sc. Returns the list and the total length.

    If binsize is given, each segment's contribution is rounded UP to the next
    bin boundary (ceil(seg_len / binsize) * binsize), so all offsets are exact
    multiples of binsize and segment boundaries coincide with bin edges.
    """
    segs = []
    offset = 0
    for seg in sc:
        if seg.assembly != assembly:
            continue
        segs.append((seg.contig, seg.start, seg.end, offset))
        seg_len = seg.end - seg.start
        if binsize is not None:
            offset += math.ceil(seg_len / binsize) * binsize
        else:
            offset += seg_len
    return segs, offset


def _bins_to_genome_intervals(bin_values: list, segs_info: list, binsize: int) -> list:
    """
    Map hot bins (bin_values[i] > 0) back to merged genome (contig, start, end) intervals.
    segs_info: [(contig, seg_start, seg_end, plot_offset), ...] from _build_axis_map.
    """
    intervals = []
    for i, v in enumerate(bin_values):
        if v <= 0:
            continue
        plot_start = i * binsize
        for contig, seg_start, seg_end, offset in segs_info:
            seg_plot_end = offset + math.ceil((seg_end - seg_start) / binsize) * binsize
            if offset <= plot_start < seg_plot_end:
                gen_start = seg_start + (plot_start - offset)
                gen_end   = min(seg_end, seg_start + ((i + 1) * binsize - offset))
                intervals.append((contig, gen_start, gen_end))
                break

    # Merge consecutive intervals on the same contig
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for contig, s, e in intervals[1:]:
        if contig == merged[-1][0] and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([contig, s, e])
    return [(c, s, e) for c, s, e in merged]


def _coord_ticks(segs_info, total_len):
    """Major and minor tick positions for genome coordinates across concatenated segments.
    Major ticks: one per segment at its start (segment boundary).
    Minor ticks: absolute genome multiples of step inside each segment; step auto-selected
    to give ≤20 ticks across total_len.
    Returns (major_positions, minor_positions) where minor excludes major positions."""
    _STEPS = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000,
              100_000, 200_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000]
    minor_step = _STEPS[-1]
    for s in _STEPS:
        if total_len / s <= 40:
            minor_step = s
            break

    major_pos = [float(off) for _, _, _, off in segs_info]
    major_set = set(major_pos)

    minor_pos = []
    for _, gen_s, gen_e, off in segs_info:
        first = math.ceil(gen_s / minor_step) * minor_step
        for g in range(first, gen_e, minor_step):
            plot_coord = off + (g - gen_s)
            if 0 < plot_coord <= total_len and float(plot_coord) not in major_set:
                minor_pos.append(float(plot_coord))

    return major_pos, minor_pos
