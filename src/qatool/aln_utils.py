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
