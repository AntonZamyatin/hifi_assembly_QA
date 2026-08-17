"""Reads-as-lines dotplot: draw read-alignment correspondences between two assemblies.

Instead of k-mer matches, each qualifying read (pid/plen above the recruiting
thresholds) contributes one line per (asm1 block, asm2 block) pair, threaded through
the read's own coordinates.  Colinear read blocks that fall on the same diagonal are
merged into a single line whose width encodes read depth.  Every block is black except
supplementary (split) blocks, which are drawn red so breakpoint-spanning evidence stands
out.  Read classification colors deliberately live only in the read-alignment plot, not here.
"""

import sqlite3

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from qatool.plotters.axis import _build_axis_map

SUPPLEMENTARY_COLOR = "#d62728"   # red: supplementary / split-read blocks

# Per-block columns needed for geometry + supplementary detection.  `pid` is last so
# callers that ignore it can keep unpacking the geometry fields positionally.
_BLOCK_COLS = ("read_name, contig, start, end, strand, "
               "soft_left, soft_right, hard_left, hard_right, q_aln_len, aln_type, pid")


def _fetch_aln_blocks(con, reads, assemblies, pid_threshold, plen_threshold, chunk=800):
    """assembly -> list of (read_name, contig, start, end, strand, sl, sr, hl, hr, q_aln_len, aln_type)."""
    out = {a: [] for a in assemblies}
    if not reads:
        return out
    ph_asm = ",".join(["?"] * len(assemblies))
    base = f"""SELECT assembly, {_BLOCK_COLS} FROM aln
               WHERE assembly IN ({ph_asm}) AND contig IS NOT NULL
                 AND start IS NOT NULL AND end IS NOT NULL
                 AND pid >= ? AND plen >= ? AND read_name IN ({{R}})"""
    for i in range(0, len(reads), chunk):
        c = reads[i:i + chunk]
        sql = base.replace("{R}", ",".join(["?"] * len(c)))
        for row in con.execute(sql, list(assemblies) + [pid_threshold, plen_threshold] + c):
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


_KEY_COLOR = {
    "supp":  SUPPLEMENTARY_COLOR,   # red: supplementary / split-read blocks
    "black": "black",               # every other block
}


def _merge_segments(segments, tol):
    """Merge colinear same-key blocks in diagonal-offset space.

    segments: list of (x0, y0, x1, y1, key).  Never merges across key, so supplementary
    (red) blocks stay separate from the backbone.  Returns list of
    (x0, y0, x1, y1, key, depth) where depth = number of blocks merged.
    """
    groups = {}
    for x0, y0, x1, y1, key in segments:
        if x1 < x0:                       # order by x so offsets are comparable
            x0, y0, x1, y1 = x1, y1, x0, y0
        forward = y1 >= y0
        offset = (y0 - x0) if forward else (y0 + x0)
        gkey = (key, forward, round(offset / tol))
        groups.setdefault(gkey, []).append((x0, x1, offset))

    out = []
    for (key, forward, _), items in groups.items():
        items.sort()
        cx0, cx1, coff, depth = items[0][0], items[0][1], items[0][2], 1
        for x0, x1, off in items[1:]:
            if x0 <= cx1:                 # overlap or abut -> extend
                cx1 = max(cx1, x1)
                coff += off               # accumulate to average later
                depth += 1
            else:
                out.append(_finish_run(cx0, cx1, coff / depth, forward, key, depth))
                cx0, cx1, coff, depth = x0, x1, off, 1
        out.append(_finish_run(cx0, cx1, coff / depth, forward, key, depth))
    return out


def _finish_run(x0, x1, offset, forward, key, depth):
    if forward:
        y0, y1 = x0 + offset, x1 + offset
    else:
        y0, y1 = offset - x0, offset - x1
    return (x0, y0, x1, y1, key, depth)


def reads_block_dotplot(sc, read_names, pid_threshold, plen_threshold, db_path,
                        ax=None, asm1="asm1", asm2="asm2", binsize=1000,
                        merge_tol=None, aln_blocks=None,
                        min_width=0.5, max_width=6.0,
                        boundary_color="black", boundary_alpha=1.0, boundary_lw=1.0):
    """Draw read-block correspondence lines (asm2 = X, asm1 = Y).

    merge_tol defaults to binsize.  aln_blocks (from _fetch_aln_blocks) is used if given,
    else the DB is queried.  Colinear blocks on the same diagonal are merged; line width
    scales with merged read depth (min_width..max_width).  Every line is black except
    supplementary blocks, which are red and never merged into the backbone.
    """
    standalone = ax is None
    if standalone:
        _, ax = plt.subplots(figsize=(16, 5))
    if merge_tol is None:
        merge_tol = binsize

    sc = sorted(sc, key=lambda s: (s.contig, s.start, s.end))
    asm1_segs, asm1_len = _build_axis_map(sc, asm1, binsize=binsize)
    asm2_segs, asm2_len = _build_axis_map(sc, asm2, binsize=binsize)
    if not asm1_segs or not asm2_segs:
        raise ValueError("SegmentCollection must contain both asm1 and asm2 segments")
    y_lookup = _build_lookup(asm1_segs)   # asm1 -> Y
    x_lookup = _build_lookup(asm2_segs)   # asm2 -> X

    reads = list(read_names)
    if aln_blocks is None:
        con = sqlite3.connect(db_path)
        try:
            aln_blocks = _fetch_aln_blocks(con, reads, [asm1, asm2], pid_threshold, plen_threshold)
        finally:
            con.close()

    # Group blocks by read, per assembly.
    by_read = {}
    for asm, blocks in aln_blocks.items():
        for b in blocks:
            by_read.setdefault(b[0], {}).setdefault(asm, []).append(b)

    segments = []   # (x0, y0, x1, y1, key)
    for rn in reads:
        per = by_read.get(rn, {})
        a1 = per.get(asm1, [])
        a2 = per.get(asm2, [])
        if not a1 or not a2:
            continue
        for (_, c1, s1, e1, st1, sl1, sr1, hl1, hr1, q1, t1, _p1) in a1:
            r0a, r1a = _read_span(st1, sl1, sr1, hl1, hr1, q1)
            v1 = _visible_read_span(y_lookup, c1, s1, e1, r0a, r1a, st1)
            if v1 is None:
                continue                  # block lands on no displayed segment at all
            for (_, c2, s2, e2, st2, sl2, sr2, hl2, hr2, q2, t2, _p2) in a2:
                r0b, r1b = _read_span(st2, sl2, sr2, hl2, hr2, q2)
                v2 = _visible_read_span(x_lookup, c2, s2, e2, r0b, r1b, st2)
                if v2 is None:
                    continue
                # Intersect the two blocks' read spans with the parts of each that are
                # actually on screen, so an overhanging block contributes its visible part.
                ro0 = max(r0a, r0b, v1[0], v2[0])
                ro1 = min(r1a, r1b, v1[1], v2[1])
                if ro1 <= ro0:
                    continue              # no shared, visible stretch of read
                ya0 = _map_in_segment(v1[2], _read_to_ref(ro0, r0a, r1a, s1, e1, st1))
                ya1 = _map_in_segment(v1[2], _read_to_ref(ro1, r0a, r1a, s1, e1, st1))
                xb0 = _map_in_segment(v2[2], _read_to_ref(ro0, r0b, r1b, s2, e2, st2))
                xb1 = _map_in_segment(v2[2], _read_to_ref(ro1, r0b, r1b, s2, e2, st2))
                supp = (t1 == "supplementary") or (t2 == "supplementary")
                segments.append((xb0, ya0, xb1, ya1, "supp" if supp else "black"))

    merged = _merge_segments(segments, merge_tol)

    # Draw black lines (supplementary in red); line width encodes merged read depth.
    if merged:
        max_depth = max(m[5] for m in merged)
        by_key = {}
        for x0, y0, x1, y1, key, depth in merged:
            lw = min_width + (max_width - min_width) * (depth / max_depth)
            by_key.setdefault(key, ([], []))
            by_key[key][0].append([(x0, y0), (x1, y1)])
            by_key[key][1].append(lw)
        # Black backbone first, red supplementary on top.
        for key in ("black", "supp"):
            if key in by_key:
                segs, widths = by_key[key]
                ax.add_collection(LineCollection(segs, colors=_KEY_COLOR[key], linewidths=widths))

    # Segment boundaries (same convention as the k-mer overlay).
    for i in range(1, len(asm2_segs)):
        same = asm2_segs[i][0] == asm2_segs[i - 1][0]
        ax.axvline(asm2_segs[i][3], color=boundary_color,
                   alpha=boundary_alpha * (0.5 if same else 1.0),
                   linewidth=boundary_lw * (0.5 if same else 1.0))
    for i in range(1, len(asm1_segs)):
        same = asm1_segs[i][0] == asm1_segs[i - 1][0]
        ax.axhline(asm1_segs[i][3], color=boundary_color,
                   alpha=boundary_alpha * (0.5 if same else 1.0),
                   linewidth=boundary_lw * (0.5 if same else 1.0))

    if standalone:
        ax.set_xlim(0, max(asm2_len, 1))
        ax.set_ylim(0, max(asm1_len, 1))
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{asm2} (concatenated segments)")
        ax.set_ylabel(f"{asm1} (concatenated segments)")
    return ax
