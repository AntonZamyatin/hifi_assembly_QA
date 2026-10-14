"""Alignment-database and BAM read queries used by the plotting functions."""

import sqlite3

import pysam

from qatool.aln_utils import cigar_stats


def _fetch_aln_from_db(con, reads: list, assemblies: list,
                       pid_threshold: float = 0.0, plen_threshold: float = 0.0,
                       chunk_size: int = 800) -> dict:
    """
    Fetch alignments from the aln table for the given reads and assemblies,
    filtering by pid >= pid_threshold and plen >= plen_threshold in SQL.
    Returns dict: assembly -> list of (read_name, contig, start, end, strand, pid, plen).
    """
    out = {asm: [] for asm in assemblies}
    if not reads:
        return out
    asm_placeholders = ",".join(["?"] * len(assemblies))
    base_sql = f"""
        SELECT assembly, read_name, contig, start, end, strand, pid, plen
        FROM aln
        WHERE assembly IN ({asm_placeholders})
          AND contig IS NOT NULL
          AND start IS NOT NULL
          AND end IS NOT NULL
          AND pid >= ?
          AND plen >= ?
          AND read_name IN ({{READS_IN}})
    """
    for i in range(0, len(reads), chunk_size):
        chunk = reads[i:i + chunk_size]
        reads_in = ",".join(["?"] * len(chunk))
        sql = base_sql.replace("{READS_IN}", reads_in)
        params = list(assemblies) + [pid_threshold, plen_threshold] + chunk
        rows = con.execute(sql, params).fetchall()
        for asm, read_name, contig, start, end, strand, pid, plen in rows:
            if pid is None or plen is None:
                raise ValueError(
                    f"Unexpected NULL pid or plen for read {read_name} in assembly {asm}"
                )
            out[str(asm)].append(
                (str(read_name), str(contig), int(start), int(end),
                 int(strand), float(pid), float(plen))
            )
    return out


def _collect_bam_read_names(bamfile: str, segments: list,
                            pid_threshold: float = 0.0,
                            plen_threshold: float | None = None) -> set:
    """Return names of all reads with qualifying alignments overlapping segments."""
    names = set()
    bam = pysam.AlignmentFile(bamfile, "rb")
    try:
        for seg in segments:
            for aln in bam.fetch(seg.contig, start=seg.start, end=seg.end):
                _, _, _, _, _, _, plen, pid, _ = cigar_stats(aln)
                if pid < pid_threshold:
                    continue
                if plen_threshold is not None and (plen is None or plen < plen_threshold):
                    continue
                names.add(aln.query_name)
    finally:
        bam.close()
    return names


# No alignment spans more reference than this (HiFi reads are tens of kb); it bounds
# _locus_read_names' index scan from below.
# ponytail: fixed bound, a read spanning more is missed; store the span if ever needed.
MAX_ALN_SPAN = 1_000_000


def _locus_read_names(con, assembly: str, segments: list,
                      pid_threshold: float, plen_threshold: float) -> set:
    """Names of reads with a qualifying alignment overlapping the segments: what
    _collect_bam_read_names finds in the BAM, from the aln_by_locus index."""
    names = set()
    for seg in segments:
        names.update(r for (r,) in con.execute(
            "SELECT read_name FROM aln WHERE assembly = ? AND contig = ? "
            "AND start > ? AND start < ? AND end > ? AND pid >= ? AND plen >= ?",
            (assembly, seg.contig, seg.start - MAX_ALN_SPAN, seg.end, seg.start,
             pid_threshold, plen_threshold)))
    return names


def _stripe_coords(db_path: str, read_set: set, assembly: str,
                   pid_threshold: float, plen_threshold: float) -> dict:
    """Return merged intervals {contig: [(start, end, n_reads), ...]} for read_set in assembly."""
    if not read_set:
        return {}
    reads_list = list(read_set)
    plen_thr = plen_threshold if plen_threshold is not None else 0.0
    intervals: dict = {}
    con = sqlite3.connect(db_path)
    try:
        for i in range(0, len(reads_list), 800):
            chunk = reads_list[i:i + 800]
            ph = ",".join(["?"] * len(chunk))
            for contig, start, end, rname in con.execute(
                f"""SELECT contig, start, end, read_name FROM aln
                    WHERE assembly = ? AND contig IS NOT NULL
                      AND start IS NOT NULL AND end IS NOT NULL
                      AND pid >= ? AND plen >= ?
                      AND read_name IN ({ph})""",
                [assembly, pid_threshold, plen_thr] + chunk,
            ):
                intervals.setdefault(str(contig), []).append((int(start), int(end), str(rname)))
    finally:
        con.close()
    merged: dict = {}
    for contig, ivs in sorted(intervals.items()):
        ivs.sort()
        cur_s, cur_e, cur_reads = ivs[0][0], ivs[0][1], {ivs[0][2]}
        out = []
        for s, e, rname in ivs[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
                cur_reads.add(rname)
            else:
                out.append((cur_s, cur_e, len(cur_reads)))
                cur_s, cur_e, cur_reads = s, e, {rname}
        out.append((cur_s, cur_e, len(cur_reads)))
        merged[contig] = out
    return merged
