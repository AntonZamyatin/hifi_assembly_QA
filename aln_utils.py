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

def calculate_pid(aln: pysam.AlignedSegment) -> float:
    """
    Estimate query percent identity from edit distance tag (NM).
    Returns:
        pid - percent of identical bases in alignment
    """
    try:
        nm = aln.get_tag("NM")
        qlen = aln.query_length
        return 1 - (nm / qlen) if qlen > 0 else 0
    except KeyError:
        return 0

def cigar_stats(aln: pysam.AlignedSegment) -> Tuple[int, int, int, int, int, int, Optional[float]]:
    """
    Calculate stats from CIGAR string.
    Returns:
      soft_left - number of softclipped bases at the alignment start
      soft_right - number of softclipped bases at the alignment end
      hard_left - number of hardclipped bases at the alignment start
      hard_right - number of hardclipped bases at the alignment end
      q_aln_len - sum of query-consuming aligned ops: M, I, =, X
      q_len - q_aln_len + soft clips + hard clips from both ends
      plen - q_aln_len / q_len  (None if q_len == 0)
    """
    ops = aln.cigartuples
    if not ops:
        return 0, 0, 0, 0, 0, 0, None

    # Left clips: allow H then S, or S then H (rare but tolerate).
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

    # Query-aligned length: query-consuming aligned ops
    # (M, I, =, X). Note: S/H are excluded by definition here.
    q_aln_len = 0
    for op, l in ops:
        if op in (0, 1, 7, 8):
            q_aln_len += l

    q_len = q_aln_len + soft_left + soft_right + hard_left + hard_right
    plen = (q_aln_len / q_len) if q_len > 0 else None
    return soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen
