"""Reads under an alignment in the right sidebar, projected base by base onto its columns.

The alignment strip (qaview.align) lays two stretches out in shared columns.
Reads aligned to the top stretch are drawn above it, reads on the bottom stretch
below, each read *projected* onto those columns through its CIGAR -- the way a
pairwise alignment extends to a multiple one:

- a base aligned to its own row goes in that row base's column: 'b' where it
  matches, 'x' where it does not;
- a row base the read skips is 'd' (a deletion relative to its row);
- where the row has gap columns -- the other assembly has bases there that this
  one lacks -- a read inserting bases at that point fills them as whole
  nucleotides, compared with the other row: 'i' where it agrees, 'x' where not;
- a read inserting nothing there has 'g' in them: it lacks those bases exactly
  as its row does;
- inserted bases with no gap columns to go in are a tick at the column after
  them. That is not "neither assembly has them": in a repeat whose length the
  assemblies disagree on, a read of a third length (say 43 bp between 42 and
  44) is one base longer than one row and one shorter than the other, and the
  lone extra base has no gap column -- those hold whole repeat units.

So at a column where the candidates differ each read shows, nucleotide for
nucleotide, which of the two it carries.

In a tandem repeat an indel can sit anywhere along the repeat at the same score,
and the aligner placed the reads' in genome orientation while the row may be
shown reverse-complemented. Each read indel is therefore walked through all its
equivalent positions in the row as shown, and put where it meets the
alignment's gap columns -- an insertion where its row has gap columns, a
deletion over bases the other row lacks -- or else at the leftmost, which is
where qaview.align left-aligns the alignment's own gaps.

The placements are exactly the ones the marginal plots count
(`region._best_placements`): pid and plen at or above the display thresholds,
and within best_pid_delta of the read's best pid in that assembly, genome-wide,
for every class -- plus the supplementary pieces of a best primary. Colour and
band are the marginal's too.

BAM records here carry no MD tag, and secondaries no SEQ, so a read's sequence
comes from its primary record (qamodel's ReadStore), is turned to the
placement's strand, and walked along the placement's CIGAR.
"""

from __future__ import annotations

import sqlite3

from qamodel.phmm.reads import ReadStore, revcomp
from qatool.aln_utils import cigar_stats, rows_for_reads

_M, _I, _D, _N, _S, _H, _EQ, _X = 0, 1, 2, 3, 4, 5, 7, 8


def column_maps(a: dict) -> tuple[list, list]:
    """Position along each row as shown -> alignment column, for top and bottom.

    Index d is the d-th base of the row left to right. For the top row and a
    forward bottom row that is genome offset d; a reverse-complemented bottom
    row runs from the stretch's end (see `_display`).
    """
    top = [k for k, c in enumerate(a["top"]) if c not in " -"]
    bottom = [k for k, c in enumerate(a["bottom"]) if c not in " -"]
    return top, bottom


def _display(g: int, lo: int, hi: int, reverse: bool) -> int:
    """Genome position -> position along the row as shown."""
    return (hi - 1 - g) if reverse else (g - lo)


def _best_info(db, assembly: str, names, pid: float, plen: float) -> dict:
    """read -> (its best non-supplementary pid, its primary's pid or None) over its
    qualifying alignments in `assembly`, genome-wide."""
    out: dict = {}
    for name, p, kind in rows_for_reads(
            db, "SELECT read_name, pid, aln_type FROM aln WHERE assembly = ? "
                "AND pid >= ? AND plen >= ? AND read_name IN ({reads})",
            names, before=[assembly, pid, plen]):
        if kind == "supplementary":
            continue
        best, prim = out.get(name, (-1.0, None))
        out[name] = (max(best, p), p if kind in ("primary", "primary_wsup") else prim)
    return out


def placements(bam, db, assembly: str, contig: str, lo: int, hi: int, *,
               pid: float, plen: float, best_pid_delta: float) -> list:
    """The alignments over [lo, hi) that the marginal for `assembly` counts: a read's
    best placements."""
    cand = []
    for a in bam.fetch(contig, lo, hi):
        if a.is_unmapped:
            continue
        stats = cigar_stats(a)
        p, pl = stats[7], stats[6]
        if p < pid or pl is None or pl < plen:
            continue
        cand.append((a, p))
    info = _best_info(db, assembly, {a.query_name for a, _ in cand}, pid, plen)
    keep = []
    for a, p in cand:
        best, prim = info.get(a.query_name, (None, None))
        if best is None:
            continue
        if a.is_supplementary:
            ok = prim is not None and prim + best_pid_delta >= best
        else:
            ok = p + best_pid_delta >= best
        if ok:
            keep.append((a, p))
    return keep


def events(a, read_seq: str, ref: str, lo: int, hi: int, reverse: bool = False) -> dict:
    """The placement's differences from its row, in row coordinates as shown.

    Returns {"span": (ds, de), "x": {d}, "d": [(d, n)], "i": [(d, bases)],
    "shown": row}: the shown positions [ds, de) the read covers, mismatches,
    deletions of n row bases from d, bases inserted before shown position d
    (strictly inside the span) as they read in the shown orientation, and the
    row's sequence as shown. Indels are where the aligner put them; `project`
    chooses among their equivalent positions.
    """
    seq = revcomp(read_seq) if a.is_reverse else read_seq
    shown = revcomp(ref) if reverse else ref
    ops = a.cigartuples
    q = ops[0][1] if ops and ops[0][0] == _H else 0      # hard-clipped bases precede
    r = a.reference_start
    s_g, e_g = max(a.reference_start, lo), min(a.reference_end, hi)
    ds, de = ((hi - e_g, hi - s_g) if reverse else (s_g - lo, e_g - lo))
    mism, dels, ins = set(), [], []
    for op, n in ops:
        if op in (_M, _EQ, _X):
            for g in range(max(r, lo), min(r + n, hi)):
                rb, qb = ref[g - lo], seq[q + (g - r)]
                if rb != qb and rb != "N" and qb != "N":
                    mism.add(_display(g, lo, hi, reverse))
            q += n
            r += n
        elif op == _I:
            if s_g < r < e_g:
                bases = seq[q:q + n]
                # before genome base r; shown reversed, before the base shown after r - 1
                d = (hi - r) if reverse else (r - lo)
                ins.append((d, revcomp(bases) if reverse else bases))
            q += n
        elif op in (_D, _N):
            s, e = max(r, lo), min(r + n, hi)
            if s < e:
                dels.append((_display(e - 1, lo, hi, True) if reverse else s - lo, e - s))
            r += n
        elif op == _S:
            q += n
    return {"span": (ds, de), "x": mism, "d": dels, "i": ins, "shown": shown}


def _place_deletion(d, n, ds, de, shown, cols, other, fixed):
    """Among the deletion's equivalent positions, one over bases the other row lacks.

    Sliding past a base is only equivalent where the read agrees with the row
    there; `fixed` holds the positions it may not cross (mismatches, other
    indels), or the alignment would change.
    """
    seen = [d]
    while d - 1 >= ds and d - 1 not in fixed and shown[d - 1] == shown[d + n - 1]:
        d -= 1
        seen.append(d)
    # also rightwards: the aligner's placement need not be the leftmost
    d = seen[0]
    while d + n < de and d + n not in fixed and shown[d] == shown[d + n]:
        d += 1
        seen.append(d)
    hits = [p for p in seen if p + n <= len(cols) and all(other[cols[k]] == "-" for k in range(p, p + n))]
    return min(hits) if hits else min(seen)


def _place_insertion(d, bases, ds, de, shown, cols, fixed):
    """Among the insertion's equivalent positions, one where its row has gap columns.

    As for deletions, it may not rotate past a position in `fixed`.
    """
    seen = [(d, bases)]
    p, b = d, bases
    while p - 1 > ds and p - 1 not in fixed and shown[p - 1] == b[-1]:   # rotate left
        b, p = b[-1] + b[:-1], p - 1
        seen.append((p, b))
    p, b = d, bases
    while p < de - 1 and p not in fixed and shown[p] == b[0]:            # and right
        b, p = b[1:] + b[0], p + 1
        seen.append((p, b))
    hits = [(p, b) for p, b in seen if cols[p] - cols[p - 1] > 1]
    return min(hits or seen)


def project(ev: dict, cols: list, other: str) -> tuple[list, list]:
    """Runs [[c0, c1, kind]] over the read's columns, and ticks [[column, n]].

    `cols` maps shown row positions to columns; `other` is the other row's
    string, for comparing bases a read inserts into this row's gap columns.
    """
    ds, de = ev["span"]
    if de <= ds:
        return [], []
    shown = ev["shown"]
    state = {}
    for d in range(ds, de):
        state[cols[d]] = "x" if d in ev["x"] else "b"
    # Positions an indel may not slide across: mismatches, and every other
    # indel's bases (deletions) or neighbours (insertions) as placed so far.
    fixed = set(ev["x"])
    for d0, n in ev["d"]:
        fixed.update(range(d0, d0 + n))
    for d, _ in ev["i"]:
        fixed.update((d - 1, d))
    placed = []
    for d0, n in ev["d"]:
        fixed.difference_update(range(d0, d0 + n))
        fixed |= ev["x"]                                  # a mismatch stays fixed
        d0 = _place_deletion(d0, n, ds, de, shown, cols, other, fixed)
        fixed.update(range(d0, d0 + n))
        placed.append((d0, n))
    for d0, n in placed:
        for d in range(max(d0, ds), min(d0 + n, de)):
            state[cols[d]] = "d"
    inserted = {}
    for d, bases in ev["i"]:
        fixed.difference_update((d - 1, d))
        fixed |= ev["x"]
        d, bases = _place_insertion(d, bases, ds, de, shown, cols, fixed)
        fixed.update((d - 1, d))
        inserted[d] = inserted.get(d, "") + bases
    ticks = []
    for d in range(ds + 1, de):
        gap = range(cols[d - 1] + 1, cols[d])             # this row's gap columns
        bases = inserted.get(d, "")
        for j, c in enumerate(gap):
            if j < len(bases):
                state[c] = "i" if bases[j] == other[c] else "x"
            else:
                state[c] = "g"
        if len(bases) > len(gap):
            ticks.append([cols[d], len(bases) - len(gap)])
    runs, c0 = [], cols[ds]
    for c in range(cols[ds] + 1, cols[de - 1] + 2):
        if c > cols[de - 1] or state[c] != state[c0]:
            runs.append([c0, c - 1, state[c0]])
            c0 = c
    return runs, ticks


def pileup(rd, a: dict, bams: dict, fastas: dict, db_path: str, stacking, *,
           pid: float, plen: float, best_pid_delta: float) -> dict:
    """Rows of reads for both sides of alignment `a` (a qaview.align result).

    `bams`/`fastas` map assembly name to a file path; `stacking` is the
    marginal's (glyphs.class_stacking or ratio_stacking). Rows are sorted by
    band, then by first column: band 0 is the marginal's bottom band.
    """
    import pysam

    top_map, bottom_map = column_maps(a)
    gidx = {n: i for i, n in enumerate(rd.read_names)}
    db = sqlite3.connect(db_path)
    bam_f = {n: pysam.AlignmentFile(p, "rb") for n, p in bams.items()}
    fa_f = {n: pysam.FastaFile(p) for n, p in fastas.items()}
    try:
        with ReadStore(db_path, bams) as store:
            sides = {}
            for side, asm, contig, lo, hi, cols, rev, other in (
                    ("top", a["top_name"], a["top_contig"], a["top_start"], a["top_end"],
                     top_map, False, a["bottom"]),
                    ("bottom", a["bottom_name"], a["bottom_contig"], a["bottom_start"],
                     a["bottom_end"], bottom_map, a["strand"] == "-", a["top"])):
                ref = fa_f[asm].fetch(contig, lo, hi).upper()
                rows = []
                for aln, p in placements(bam_f[asm], db, asm, contig, lo, hi, pid=pid, plen=plen,
                                         best_pid_delta=best_pid_delta):
                    g = gidx.get(aln.query_name)
                    band = stacking.read_band[g] if g is not None else len(stacking.keys) - 1
                    try:
                        seq, _ = store.get(aln.query_name)
                    except (KeyError, ValueError):
                        seq = None                         # no stored sequence
                    if seq is not None:
                        runs, ticks = project(events(aln, seq, ref, lo, hi, rev), cols, other)
                    else:                                  # the span alone, as one bar
                        s, e = max(aln.reference_start, lo), min(aln.reference_end, hi)
                        c0, c1 = sorted((cols[_display(s, lo, hi, rev)],
                                         cols[_display(e - 1, lo, hi, rev)]))
                        runs, ticks = [[c0, c1, "b"]], []
                    if not runs:
                        continue
                    rows.append({"name": aln.query_name, "band": int(band),
                                 "col0": runs[0][0], "col1": runs[-1][1], "pid": round(p, 5),
                                 "start": aln.reference_start, "end": aln.reference_end,
                                 "strand": "-" if aln.is_reverse else "+",
                                 "supplementary": aln.is_supplementary,
                                 "sequence": seq is not None, "runs": runs, "ticks": ticks})
                rows.sort(key=lambda row: (row["band"], row["col0"], row["col1"]))
                sides[side] = rows
    finally:
        db.close()
        for f in (*bam_f.values(), *fa_f.values()):
            f.close()
    return {"bands": [{"name": k, "color": c} for k, c in zip(stacking.keys, stacking.colors)],
            **sides}
