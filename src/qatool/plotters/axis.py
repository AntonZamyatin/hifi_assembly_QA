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


def contig_runs(segs_info, total):
    """Yield (contig, bar_off, bar_end) merging consecutive same-contig segments.

    The recruiter groups same-contig segments together, so a run is the span a single
    contig colour bar should cover.
    """
    if not segs_info:
        return
    cur_contig, _, _, cur_off = segs_info[0]
    for i in range(1, len(segs_info)):
        if segs_info[i][0] != cur_contig:
            yield cur_contig, cur_off, segs_info[i][3]
            cur_contig, _, _, cur_off = segs_info[i]
    yield cur_contig, cur_off, total


def plot_to_genome(segs_info, pos):
    """
    Plot coordinate -> (contig, genomic position); inverse of the offset mapping
    built by _build_axis_map.

    Returns None when pos falls outside every segment.  With a binsize-aligned axis
    map that includes the padding after a segment whose length is not a whole number
    of bins: those coordinates have no genomic counterpart.
    """
    for contig, start, end, offset in segs_info:
        if offset <= pos < offset + (end - start):
            return contig, start + int(pos - offset)
    return None


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
