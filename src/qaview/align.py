"""Align the stretches under a rectangle drawn on the dotplot, for the right sidebar.

The vertical axis is asm1 and the horizontal asm2. Shift+drag draws a rectangle;
each side is snapped outward to whole display bins, and the two stretches are
shown whole, laid out along their true correspondence, so any difference
between the candidates inside them is visible.

The stretches are usually offset, and aligning them end to end does not show
the offset as end gaps: under edit distance a mismatch is cheaper than two
gaps, so offset stretches collapse into a noisy alignment that looks like
divergence. Instead:

1. The asm1 stretch is located on asm2 (edlib's HW mode: free end gaps in the
   target), in both orientations, within a window centred on the rectangle: if
   the correspondence passes through the rectangle's centre, the counterpart
   of an asm1 stretch of length L lies within L of the centre along asm2, so
   that span plus half a bin is searched, together with the asm2 stretch
   itself. Centring it on the asm2 bins instead missed the counterpart when
   the rectangle sat near a bin edge, and in a repeat the fit then landed on a
   neighbouring copy.
2. The path edlib returns is polished: under its unit costs one 12 bp indel and
   six scattered small gaps with chance matches between them cost the same.
   Every stretch holding gaps is realigned locally with affine gap costs
   (parasail NW), so one event is drawn as one gap. Doing it locally keeps
   memory linear, where a full affine table for 200 kb a side would not fit.
3. Every gap is shifted as far left as the sequence allows, in the orientation
   shown. In a tandem repeat an indel can sit anywhere along the repeat at the
   same score; fixing it at the leftmost place on screen lets the pileup put
   the reads' own indels (normalised the same way, qaview.pileup) in the same
   column as the alignment's.
4. The rows show the whole asm1 stretch and the whole asm2 stretch along that
   path. Where both have sequence is the comparison; sequence of one whose
   counterpart lies outside the other is shown against a blank ("outside") and
   is not counted as a difference.

When the reverse orientation fits better it is the asm2 side that is
reverse-complemented, so the top row always reads asm1 forward.
"""

from __future__ import annotations

import edlib
import parasail

from qatool.aln_utils import revcomp

MAX_SPAN = 200_000                  # bp per side

# Affine scoring for the displayed alignment: +2 match, -4 mismatch, a gap of
# length k costs 6 + (k - 1). N matches nothing.
_MATRIX = parasail.matrix_create("ACGTN", 2, -4)
_GAP_OPEN, _GAP_EXTEND = 6, 1
_PAD = 20                           # matching columns kept around a stretch being polished


def span_at(segs_info, binsize: int, p0: float, p1: float, flips=None):
    """The bin-snapped genome stretch under plot interval [p0, p1] on one axis.

    Returns (contig, seg_start, seg_end, lo, hi, gmid, clipped). Bins restart at
    every segment's genome start, and a stretch never leaves one segment: it is
    the segment under the interval's midpoint, gmid is that midpoint in genome
    coordinates, and `clipped` says whether the interval reached beyond the
    segment. `flips` flags segments drawn backwards (`region.arrange`).
    """
    p0, p1 = sorted((p0, p1))
    mid = (p0 + p1) / 2
    for i, (contig, start, end, offset) in enumerate(segs_info):
        length = end - start
        if mid < offset + length or i == len(segs_info) - 1:
            break
    clipped = p0 < offset or p1 > offset + length
    if flips and flips[i]:                     # axis runs end -> start here
        p0, p1, mid = (2 * offset + length - p for p in (p1, p0, mid))
    rel0 = int(min(max(p0 - offset, 0), length - 1))
    rel1 = int(min(max(p1 - offset, 0), length - 1))
    lo = start + (rel0 // binsize) * binsize
    hi = min(end, start + (rel1 // binsize + 1) * binsize)
    gmid = start + int(min(max(mid - offset, 0), length - 1))
    return contig, start, end, lo, hi, gmid, clipped


def _expand(cigar: str) -> list:
    out, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            out.extend(ch * int(num))
            num = ""
    return out


def _affine(query: str, target: str) -> list:
    if not query:
        return ["D"] * len(target)
    if not target:
        return ["I"] * len(query)
    r = parasail.nw_trace_scan_16(query, target, _GAP_OPEN, _GAP_EXTEND, _MATRIX)
    if r.saturated:
        r = parasail.nw_trace_scan_32(query, target, _GAP_OPEN, _GAP_EXTEND, _MATRIX)
    return _expand(r.cigar.decode.decode())


def polish(ops: list, query: str, target: str) -> list:
    """Realign every stretch of `ops` holding a gap with affine costs; one op per column.

    Non-match columns closer than 2 * _PAD are one stretch; each stretch is
    realigned with _PAD matching columns on either side, and spliced back.
    """
    bad = [k for k, op in enumerate(ops) if op != "="]
    if not bad:
        return ops
    stretches, s0, s1 = [], bad[0], bad[0]
    for k in bad[1:]:
        if k - s1 <= 2 * _PAD:
            s1 = k
        else:
            stretches.append((s0, s1))
            s0 = s1 = k
    stretches.append((s0, s1))

    # query/target position before each column
    qpos, tpos, i, j = [], [], 0, 0
    for op in ops:
        qpos.append(i)
        tpos.append(j)
        i += op in "=XI"
        j += op in "=XD"
    qpos.append(i)
    tpos.append(j)

    out, last = [], 0
    for s0, s1 in stretches:
        if not any(ops[k] in "ID" for k in range(s0, s1 + 1)):
            continue                                      # mismatches only: nothing to move
        c0, c1 = max(0, s0 - _PAD), min(len(ops), s1 + 1 + _PAD)
        out.extend(ops[last:c0])
        out.extend(_affine(query[qpos[c0]:qpos[c1]], target[tpos[c0]:tpos[c1]]))
        last = c1
    out.extend(ops[last:])
    return out


def layout(ops: list, query: str, tgt: str, a0: int, b0: int, b1: int) -> dict:
    """Rows for the whole query and the whole stretch tgt[b0:b1], along a path.

    `ops` aligns all of `query` to tgt from a0 on, one op per column. Target
    bases outside [b0, b1) are not part of the stretch and are left out; stretch
    bases the path does not reach are added at the ends. Columns where one row
    has no stretch sequence carry a blank there -- '-' is kept for real gaps.
    """
    top, mid, bot = [], [], []
    counts = {"matches": 0, "mismatches": 0, "ins": 0, "ins_opens": 0,
              "del": 0, "del_opens": 0, "outside_top": 0, "outside_bottom": 0}

    def col(t, m, b):
        top.append(t)
        mid.append(m)
        bot.append(b)

    for k in range(b0, min(a0, b1)):                      # stretch bases before the path
        col(" ", " ", tgt[k])
        counts["outside_bottom"] += 1
    i, j, prev = 0, a0, None
    for op in ops:
        if op in "=X":
            if b0 <= j < b1:
                col(query[i], "|" if op == "=" else ".", tgt[j])
                counts["matches" if op == "=" else "mismatches"] += 1
            else:
                col(query[i], " ", " ")
                counts["outside_top"] += 1
            i += 1
            j += 1
        elif op == "I":                                   # a base of the query only
            if b0 < j < b1:
                col(query[i], " ", "-")
                counts["ins"] += 1
                counts["ins_opens"] += prev != "I"
            else:
                col(query[i], " ", " ")
                counts["outside_top"] += 1
            i += 1
        elif op == "D":                                   # a base of the target only
            if b0 <= j < b1:
                col("-", " ", tgt[j])
                counts["del"] += 1
                counts["del_opens"] += prev != "D"
            j += 1                                        # outside the stretch: not shown
        else:
            raise ValueError(f"unexpected alignment op {op!r}")
        prev = op
    for k in range(max(j, b0), b1):                       # stretch bases after the path
        col(" ", " ", tgt[k])
        counts["outside_bottom"] += 1

    compared = counts["matches"] + counts["mismatches"] + counts["ins"] + counts["del"]
    counts["columns"] = len(mid)
    counts["compared"] = compared
    counts["identity"] = counts["matches"] / compared if compared else 0.0
    return {"top": "".join(top), "mid": "".join(mid), "bottom": "".join(bot), **counts}


def left_align(top: str, mid: str, bot: str) -> tuple[str, str, str]:
    """Shift every gap as far left as the flanking sequence allows.

    A gap run in one row moves left one column while the base before it in that
    row matches the other row's base both there and at the run's last column, so
    each shift swaps one match for another and the alignment's score and counts
    are unchanged. Blank ("outside") columns stop it.
    """
    t, m, b = list(top), list(mid), list(bot)
    for gap, other in ((t, b), (b, t)):
        k = 0
        while k < len(gap):
            if gap[k] != "-":
                k += 1
                continue
            g0 = k
            while k < len(gap) and gap[k] == "-":
                k += 1
            g1 = k                                        # the run is [g0, g1)
            while (g0 > 0 and m[g0 - 1] == "|" and gap[g0 - 1] == other[g1 - 1]
                   and other[g1 - 1] not in " -"):
                gap[g1 - 1], gap[g0 - 1] = gap[g0 - 1], "-"
                m[g1 - 1], m[g0 - 1] = "|", " "
                g0 -= 1
                g1 -= 1
    return "".join(t), "".join(m), "".join(b)


def _runs(row: str, mid: str, ch: str) -> int:
    """Gap runs of `ch` in `row` at columns where the other row has a base."""
    n, inside = 0, False
    for c, x in zip(row, mid):
        on = c == ch
        n += on and not inside
        inside = on
    return n


def align_rect(rd, fasta1, fasta2, x0: float, x1: float, y0: float, y1: float) -> dict:
    """The alignment for a rectangle [x0, x1] x [y0, y1] in plot coordinates.

    A degenerate rectangle (a shift+click) aligns the single bins under it.
    `fasta1`/`fasta2` need `.fetch`. Raises ValueError for a side over MAX_SPAN.
    """
    b = rd.binsize
    flips = getattr(rd, "flips", None) or {}
    c1, _, _, lo1, hi1, _, clip1 = span_at(rd.asm1_segs_info, b, y0, y1, flips.get(rd.asm1))
    c2, s2, e2, lo2, hi2, gx, clip2 = span_at(rd.asm2_segs_info, b, x0, x1, flips.get(rd.asm2))
    for name, lo, hi in ((rd.asm1, lo1, hi1), (rd.asm2, lo2, hi2)):
        if hi - lo > MAX_SPAN:
            raise ValueError(f"the {name} side covers {hi - lo:,} bp; "
                             f"at most {MAX_SPAN:,} bp a side can be aligned")
    reach = (hi1 - lo1) + b // 2
    t_lo, t_hi = max(s2, min(lo2, gx - reach)), min(e2, max(hi2, gx + reach))
    query = fasta1.fetch(c1, lo1, hi1).upper()
    window = fasta2.fetch(c2, t_lo, t_hi).upper()

    # Pick the orientation on distance alone, then trace only the winner.
    fits = [(edlib.align(query, tgt, mode="HW", task="distance")["editDistance"], strand, tgt)
            for strand, tgt in (("+", window), ("-", revcomp(window)))]
    _, strand, tgt = min(f for f in fits if f[0] >= 0)
    r = edlib.align(query, tgt, mode="HW", task="path")
    a0, a1 = r["locations"][0]                            # inclusive, in tgt
    ops = polish(_expand(r["cigar"]), query, tgt[a0:a1 + 1])
    if strand == "+":                                     # the asm2 stretch, in tgt coordinates
        b0, b1 = lo2 - t_lo, hi2 - t_lo
        f0, f1 = t_lo + a0, t_lo + a1 + 1
    else:                                                 # tgt runs t_hi -> t_lo
        b0, b1 = t_hi - hi2, t_hi - lo2
        f0, f1 = t_hi - a1 - 1, t_hi - a0
    out = layout(ops, query, tgt, a0, b0, b1)
    out["top"], out["mid"], out["bottom"] = left_align(out["top"], out["mid"], out["bottom"])
    out["ins_opens"] = _runs(out["bottom"], out["mid"], "-")   # top-only bases: gaps below
    out["del_opens"] = _runs(out["top"], out["mid"], "-")
    out.update(
        strand=strand, binsize=b, clipped=[clip1, clip2],
        top_name=rd.asm1, top_contig=c1, top_start=lo1, top_end=hi1,
        bottom_name=rd.asm2, bottom_contig=c2, bottom_start=lo2, bottom_end=hi2,
        counterpart=[f0, f1],               # where the asm1 stretch lands on asm2
    )
    return out
