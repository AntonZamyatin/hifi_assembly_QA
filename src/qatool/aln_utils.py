"""
Project: hifi assembly QA
Module: aln_utils.py
Description: Utilities for dealing with alignments using pysam

Author: Anton Zamyatin
Affiliation: Pennsylvania State University
Contact: zamyatin.anton@gmail.com

Created: 2/11/26

License: MIT License
"""

import pysam
from typing import Tuple, Optional


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def same_sequences(loci1: list, loci2: list, seq1, seq2) -> bool:
    """The common-read test: True iff the sequences at loci1 (in one assembly) and
    loci2 (in the other) are the same multiset, each taken up to reverse
    complement. loci* = [(contig, start, end), ...]; seq*(contig, start, end)
    returns the upper-case sequence there."""
    if len(loci1) != len(loci2) or \
            sorted(e - s for _, s, e in loci1) != sorted(e - s for _, s, e in loci2):
        return False
    if len(loci1) == 1:
        a, b = seq1(*loci1[0]), seq2(*loci2[0])
        return a == b or a == revcomp(b)

    def canon(seq, loci):
        return sorted(min(x, revcomp(x)) for x in (seq(*locus) for locus in loci))
    return canon(seq1, loci1) == canon(seq2, loci2)


def rows_for_reads(con, sql: str, reads, before=(), after=(), chunk: int = 800):
    """Rows of `sql` over every read in `reads`. Its `{reads}` stands for the
    read_name IN list, run in chunks under SQLite's limit on bound variables;
    `before` and `after` bind the placeholders on either side of it."""
    reads = list(reads)
    for i in range(0, len(reads), chunk):
        c = reads[i:i + chunk]
        yield from con.execute(sql.replace("{reads}", ",".join("?" * len(c))), [*before, *c, *after])

def cigar_stats(aln: pysam.AlignedSegment) -> Tuple[int, int, int, int, int, int, Optional[float], float, int]:
    """
    Calculate stats from CIGAR string.
    Returns:
      soft_left  - softclipped bases at the alignment start
      soft_right - softclipped bases at the alignment end
      hard_left  - hardclipped bases at the alignment start
      hard_right - hardclipped bases at the alignment end
      q_aln_len  - query-consuming aligned bases: M + I + = + X
      q_len      - q_aln_len + all clips
      plen       - q_aln_len / q_len  (None if q_len == 0)
      pid        - percent identity [0,1]:
                   extended CIGAR (=,X present): = / (= + X + I + D)
                   regular CIGAR (M only):       1 - NM / q_aln_len
      ed         - edit distance:
                   extended CIGAR: X + I + D (computed from ops)
                   regular CIGAR:  NM tag value
    """
    ops = aln.cigartuples
    if not ops:
        return 0, 0, 0, 0, 0, 0, None, 0.0, 0

    # Left clips: H then S, S then H, or either alone.
    i = 0
    soft_left = hard_left = 0
    while i < len(ops) and ops[i][0] in (4, 5):
        op, l = ops[i]
        if op == 4:
            soft_left += l
        else:
            hard_left += l
        i += 1

    # Right clips
    j = len(ops) - 1
    soft_right = hard_right = 0
    while j >= 0 and ops[j][0] in (4, 5):
        op, l = ops[j]
        if op == 4:
            soft_right += l
        else:
            hard_right += l
        j -= 1

    # Accumulate per-op counts in one pass.
    matches = mismatches = insertions = deletions = m_bases = 0
    for op, l in ops:
        if op == 7:   # = sequence match
            matches += l
        elif op == 8: # X sequence mismatch
            mismatches += l
        elif op == 1: # I insertion to reference
            insertions += l
        elif op == 2: # D deletion from reference
            deletions += l
        elif op == 0: # M match-or-mismatch (no extended CIGAR)
            m_bases += l

    q_aln_len = matches + mismatches + insertions + m_bases  # query-consuming aligned ops
    q_len = q_aln_len + soft_left + soft_right + hard_left + hard_right
    plen = (q_aln_len / q_len) if q_len > 0 else None

    if matches + mismatches > 0:
        # Extended CIGAR (= and X present): compute directly from ops.
        ed = mismatches + insertions + deletions
        aln_len = matches + mismatches + insertions + deletions
        pid = matches / aln_len if aln_len > 0 else 0.0
    else:
        # Regular M CIGAR: fall back to NM tag.
        nm = int(aln.get_tag("NM")) if aln.has_tag("NM") else 0
        ed = nm
        pid = (1.0 - nm / q_aln_len) if q_aln_len > 0 else 0.0

    return soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid, ed
