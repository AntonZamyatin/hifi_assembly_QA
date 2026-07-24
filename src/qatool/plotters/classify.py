"""Read classification: support categories (asm1/asm2/common/…) and occurrence buckets."""

import hashlib
import sqlite3

import pysam

from qatool.plotters.colors import BEST_PID_DELTA_DEFAULT, PID_TIE_EPS_DEFAULT
from qatool.plotters.db import _fetch_aln_from_db


def _interval_in_segments(seg_lookup_by_contig: dict, contig: str, start: int, end: int) -> bool:
    """seg_lookup_by_contig: {contig: [(start, end), ...]} — pre-built from segs_by_asm."""
    for s, e in seg_lookup_by_contig.get(contig, ()):
        if start < e and end > s:
            return True
    return False


def _is_read_common(loci1: list, loci2: list,
                    asm1_fasta: pysam.FastaFile, asm2_fasta: pysam.FastaFile) -> bool:
    """
    Mirror of Recruiter._read_is_common: returns True iff the reference sequences
    at loci1 (asm1) and loci2 (asm2) are identical sets (including reverse complements).
    loci* = [(contig, start, end), ...]
    """
    if len(loci1) != len(loci2):
        return False
    _comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    def _hashes(fa, loci):
        c: dict = {}
        for contig, s, e in loci:
            seq = fa.fetch(contig, s, e).upper()
            h1 = hashlib.sha1(seq.encode()).digest()
            h2 = hashlib.sha1(seq.translate(_comp)[::-1].encode()).digest()
            c[h1] = c.get(h1, 0) + 1
            c[h2] = c.get(h2, 0) + 1
        return c
    return _hashes(asm1_fasta, loci1) == _hashes(asm2_fasta, loci2)


def _classify_reads_for_occurrence_colors(
    aln_by_asm: dict,
    asm: str,
    common_reads: set,
    best_pid_delta: float = BEST_PID_DELTA_DEFAULT,
) -> list:
    """
    Bucket reads by number of best-PID alignments to asm.

    Returns 10 sets matching OCCURRENCE_COLORS:
        indices 0-4 : 1 / 2 / 3 / 4 / 5+ best alignments, non-common reads
        indices 5-9 : same counts for common reads (pale colors)
    """
    alns = aln_by_asm.get(asm, [])

    best_pid: dict = {}
    for read_name, contig, start, end, strand, pid, plen in alns:
        if pid > best_pid.get(read_name, float("-inf")):
            best_pid[read_name] = pid

    n_best: dict = {}
    for read_name, contig, start, end, strand, pid, plen in alns:
        if pid + best_pid_delta >= best_pid.get(read_name, float("-inf")):
            n_best[read_name] = n_best.get(read_name, 0) + 1

    buckets: list = [set() for _ in range(10)]
    for read_name, count in n_best.items():
        idx = min(count - 1, 4) * 2          # 0, 2, 4, 6, 8
        if read_name in common_reads:
            idx += 1                           # 1, 3, 5, 7, 9
        buckets[idx].add(read_name)

    return buckets


def _classify_reads_for_support_colors(db_file: str, reads_set: set[str],
                                       pid_thr: float, plen_thr: float | None,
                                       segs_by_asm: dict | None, asm1: str, asm2: str,
                                       pid_tie_eps: float = PID_TIE_EPS_DEFAULT,
                                       aln_by_asm: dict | None = None,
                                       asm1_fasta: pysam.FastaFile | None = None,
                                       asm2_fasta: pysam.FastaFile | None = None,
                                       best_pid_delta: float = BEST_PID_DELTA_DEFAULT) -> dict:
    """
    Classify each read into: asm1_only, asm2_only, asm1_best, asm2_best, common, other.

    For reads aligning to both assemblies within the recruited segments:
      - common  (gold): reference sequences at both loci are identical (SHA1 including RC),
                        matching Recruiter._read_is_common.  Requires asm1_fasta/asm2_fasta;
                        without them any equal-pid dual-mapper is labelled common.
      - asm1_best / asm2_best: pid difference exceeds eps (and sequences are not identical).
      - other: pid tie but sequences differ (or no qualifying aln in either assembly).

    If aln_by_asm is provided it is used directly and no DB query is issued.
    """
    reads = list(reads_set)
    out = {
        "asm1_only": set(),
        "asm2_only": set(),
        "asm1_best": set(),
        "asm2_best": set(),
        "common":    set(),
        "other":     set(),
    }
    if not reads:
        return out

    if aln_by_asm is None:
        con = sqlite3.connect(db_file)
        try:
            aln_by_asm = _fetch_aln_from_db(
                con, reads, [asm1, asm2],
                pid_threshold=pid_thr,
                plen_threshold=plen_thr if plen_thr is not None else 0.0,
            )
        finally:
            con.close()

    seg_lookup: dict | None = None
    if segs_by_asm is not None:
        seg_lookup = {}
        for asm, seg_list in segs_by_asm.items():
            by_contig: dict = {}
            for seg in seg_list:
                by_contig.setdefault(seg.contig, []).append((seg.start, seg.end))
            seg_lookup[asm] = by_contig

    # Pass 1: best PID per (read, assembly), respecting the segment filter.
    best_pid: dict = {}
    for asm, alns in aln_by_asm.items():
        for read_name, contig, start, end, strand, pid, plen in alns:
            if seg_lookup is not None and not _interval_in_segments(seg_lookup.get(asm, {}), contig, start, end):
                continue
            key = (read_name, asm)
            if pid > best_pid.get(key, float("-inf")):
                best_pid[key] = pid

    # Pass 2: admit only alignments within best_pid_delta of the per-(read,assembly) best.
    per_read: dict = {}
    for asm, alns in aln_by_asm.items():
        for read_name, contig, start, end, strand, pid, plen in alns:
            if seg_lookup is not None and not _interval_in_segments(seg_lookup.get(asm, {}), contig, start, end):
                continue
            if pid + best_pid_delta < best_pid.get((read_name, asm), float("-inf")):
                continue
            info = per_read.setdefault(read_name, {}).setdefault(
                asm, {"n_aln": 0, "max_pid": float("-inf"), "loci": []}
            )
            info["n_aln"] += 1
            info["max_pid"] = max(info["max_pid"], pid)
            info["loci"].append((contig, start, end))

    for rn in reads:
        asm_vals = per_read.get(rn, {})
        asm1_info = asm_vals.get(asm1)
        asm2_info = asm_vals.get(asm2)
        has_asm1 = asm1_info is not None and asm1_info["n_aln"] > 0
        has_asm2 = asm2_info is not None and asm2_info["n_aln"] > 0
        if has_asm1 and not has_asm2:
            out["asm1_only"].add(rn)
        elif has_asm2 and not has_asm1:
            out["asm2_only"].add(rn)
        elif has_asm1 and has_asm2:
            d = asm2_info["max_pid"] - asm1_info["max_pid"]
            pids_tied = abs(d) <= best_pid_delta
            if pids_tied and asm1_fasta is not None and asm2_fasta is not None:
                if _is_read_common(asm1_info["loci"], asm2_info["loci"], asm1_fasta, asm2_fasta):
                    out["common"].add(rn)
                else:
                    if d > pid_tie_eps:
                        out["asm2_best"].add(rn)
                    elif -d > pid_tie_eps:
                        out["asm1_best"].add(rn)
                    # else: different loci, equal PID → "other"
            elif pids_tied:
                # no fastas but PID tie → assume common
                out["common"].add(rn)
            else:
                # PIDs differ by more than best_pid_delta → classify by PID, skip hash check
                if d > pid_tie_eps:
                    out["asm2_best"].add(rn)
                elif -d > pid_tie_eps:
                    out["asm1_best"].add(rn)
                # else: abs(d) > best_pid_delta but within pid_tie_eps → "other"

    classified = out["asm1_only"] | out["asm2_only"] | out["asm1_best"] | out["asm2_best"] | out["common"]
    out["other"] = set(reads_set) - classified
    return out
