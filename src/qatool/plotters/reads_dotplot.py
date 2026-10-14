"""Read-alignment block geometry: the aln rows qaview draws as read lines, and the
maps between read, reference and plot coordinates (qaview.region, reads_overlay).
"""

from qatool.aln_utils import rows_for_reads

# Per-block columns needed for geometry + supplementary detection.  `pid` is last so
# callers that ignore it can keep unpacking the geometry fields positionally.
_BLOCK_COLS = ("read_name, contig, start, end, strand, "
               "soft_left, soft_right, hard_left, hard_right, q_aln_len, aln_type, pid")


def _fetch_aln_blocks(con, reads, assemblies, pid_threshold, plen_threshold):
    """assembly -> list of (read_name, contig, start, end, strand, sl, sr, hl, hr, q_aln_len, aln_type)."""
    out = {a: [] for a in assemblies}
    if not reads:
        return out
    ph_asm = ",".join(["?"] * len(assemblies))
    base = f"""SELECT assembly, {_BLOCK_COLS} FROM aln
               WHERE assembly IN ({ph_asm}) AND contig IS NOT NULL
                 AND start IS NOT NULL AND end IS NOT NULL
                 AND pid >= ? AND plen >= ? AND read_name IN ({{reads}})"""
    for row in rows_for_reads(con, base, reads, before=[*assemblies, pid_threshold, plen_threshold]):
        out[str(row[0])].append(row[1:])
    return out


def _read_span(strand, sl, sr, hl, hr, q_aln_len):
    """Aligned read portion in original (5'->3') read coordinates.

    For a reverse-strand record the clips reported on the alignment's left are the
    read's 3' side, so the 5' offset comes from the right clips instead.
    """
    r0 = (hl + sl) if strand == 1 else (hr + sr)
    return r0, r0 + q_aln_len


def _read_to_ref(rp, r0, r1, ref_start, ref_end, strand):
    """Linear read-position -> reference-position within one block."""
    if r1 <= r0:
        return float(ref_start)
    frac = (rp - r0) / (r1 - r0)
    if strand == 1:
        return ref_start + frac * (ref_end - ref_start)
    return ref_end - frac * (ref_end - ref_start)          # reverse: 5' read -> ref_end


def _build_lookup(segs):
    """contig -> [(seg_start, seg_end, plot_offset), ...] from _build_axis_map output."""
    lk = {}
    for contig, s, e, off in segs:
        lk.setdefault(contig, []).append((s, e, off))
    return lk


def _map_point(lookup, contig, ref_start, ref_end, pos):
    """Map a reference position to plot coord, if [ref_start,ref_end] fits one displayed segment.

    Strict containment.  Kept for callers that need the old all-or-nothing test;
    the drawing paths use `_visible_read_span` + `_map_in_segment` instead, which clip.
    """
    for s, e, off in lookup.get(contig, ()):
        if s <= ref_start and ref_end <= e:
            return off + (pos - s)
    return None


def _ref_to_read(ref_pos, r0, r1, ref_start, ref_end, strand):
    """Inverse of `_read_to_ref` within one block."""
    span = ref_end - ref_start
    if span <= 0:
        return r0
    frac = (ref_pos - ref_start) / span if strand == 1 else (ref_end - ref_pos) / span
    return r0 + frac * (r1 - r0)


def _overlapping_segment(lookup, contig, ref_start, ref_end):
    """Displayed segment overlapping [ref_start, ref_end] most, or None.

    Replaces a full-containment test: requiring a block to sit entirely inside one
    segment discards it whole the moment it pokes past an edge, which loses
    2L/(S+L) of reads for read length L and segment length S — about half at 15 kb
    reads over 41 kb segments, and everything once segments drop below read length.
    """
    best, best_ov = None, 0
    for s, e, off in lookup.get(contig, ()):
        ov = min(ref_end, e) - max(ref_start, s)
        if ov > best_ov:
            best, best_ov = (s, e, off), ov
    return best


def _visible_read_span(lookup, contig, ref_start, ref_end, r0, r1, strand):
    """(lo, hi, segment): the read-coordinate interval whose reference image is on screen.

    Clipping in read space (rather than clipping each axis independently) is what
    keeps a line's slope honest: both ends still refer to the same stretch of read.
    """
    seg = _overlapping_segment(lookup, contig, ref_start, ref_end)
    if seg is None:
        return None
    s, e, _ = seg
    vs, ve = max(ref_start, s), min(ref_end, e)
    if ve <= vs:
        return None
    a = _ref_to_read(vs, r0, r1, ref_start, ref_end, strand)
    b = _ref_to_read(ve, r0, r1, ref_start, ref_end, strand)
    return (a, b, seg) if a <= b else (b, a, seg)


def _map_in_segment(seg, ref_pos):
    """Reference position -> plot coordinate within an already-chosen segment."""
    s, e, off = seg
    return off + min(max(ref_pos, s), e) - s
