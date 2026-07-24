"""
Project: hifi assembly QA
Module: analysis.py

Author: Anton Zamyatin
Affiliation: Pennsylvania State University
Contact: zamyatin.anton@gmail.com
"""

import re
import bisect
import sqlite3
import pysam
import pandas as pd
from typing import Dict, List, Optional, Tuple


def parse_stripe_coords(s: str) -> Dict[str, List[Tuple[int, int, int]]]:
    """
    Parse stripe coordinate string from likelihood_report output into a dict.

    Input:  "[ctg1:100000-105000(8), ctg1:210000-211500(3)]"
    Output: {"ctg1": [(100000, 105000, 8), (210000, 211500, 3)]}
    """
    result: Dict[str, List[Tuple[int, int, int]]] = {}
    for m in re.finditer(r'(\S+?):(\d+)-(\d+)\((\d+)\)', s):
        contig = m.group(1)
        result.setdefault(contig, []).append(
            (int(m.group(2)), int(m.group(3)), int(m.group(4)))
        )
    return result



def reads_aln_df(
    db_path: str,
    read_names,
    assembly: str,
) -> pd.DataFrame:
    """
    Return all alignments (no quality thresholds) to assembly for the given reads.

    Parameters
    ----------
    db_path : str
    read_names : iterable of str  (e.g. df["read_name"].unique())
    assembly : str

    Returns
    -------
    pd.DataFrame — same columns as blue_stripe_aln_df.
    Unmapped rows (contig IS NULL) are included.
    """
    reads_list = list(read_names)
    if not reads_list:
        return pd.DataFrame()

    rows_out = []
    con = sqlite3.connect(db_path)
    try:
        for i in range(0, len(reads_list), 800):
            chunk = reads_list[i:i + 800]
            ph = ",".join(["?"] * len(chunk))
            rows_out.extend(con.execute(
                f"""SELECT read_name, contig, start, end, strand,
                           pid, plen, ed, q_aln_len, q_len,
                           soft_left, soft_right, aln_type
                    FROM aln
                    WHERE assembly = ?
                      AND read_name IN ({ph})""",
                [assembly] + chunk,
            ).fetchall())
    finally:
        con.close()

    df = pd.DataFrame(rows_out, columns=[
        "read_name", "contig", "start", "end", "strand",
        "pid", "plen", "ed", "q_aln_len", "q_len",
        "soft_left", "soft_right", "aln_type",
    ])
    return df.sort_values(["contig", "start"]).reset_index(drop=True)


# ── region / read queries ─────────────────────────────────────────────────────

_ALN_COLS = [
    "read_name", "contig", "start", "end", "strand",
    "pid", "plen", "ed", "q_aln_len", "q_len",
    "soft_left", "soft_right", "aln_type",
]


def segments_aln_df(
    db_path: str,
    sc,
    assembly: str,
) -> pd.DataFrame:
    """
    Return all alignments that overlap any segment in sc for the given assembly.
    No pid/plen filtering — every row in the DB that overlaps is returned.

    Parameters
    ----------
    db_path : str
    sc : SegmentCollection (iterable of Segment)
        Only segments whose .assembly == assembly are queried.
    assembly : str

    Returns
    -------
    pd.DataFrame with columns: read_name, contig, start, end, strand,
        pid, plen, ed, q_aln_len, q_len, soft_left, soft_right, aln_type.
    Sorted by contig, start.
    """
    segs = [s for s in sc if s.assembly == assembly]
    if not segs:
        return pd.DataFrame(columns=_ALN_COLS)

    rows_out = []
    con = sqlite3.connect(db_path)
    try:
        for seg in segs:
            rows_out.extend(con.execute(
                """SELECT read_name, contig, start, end, strand,
                          pid, plen, ed, q_aln_len, q_len,
                          soft_left, soft_right, aln_type
                   FROM aln
                   WHERE assembly = ? AND contig = ?
                     AND start < ? AND end > ?""",
                [assembly, seg.contig, seg.end, seg.start],
            ).fetchall())
    finally:
        con.close()

    df = pd.DataFrame(rows_out, columns=_ALN_COLS)
    return df.drop_duplicates().sort_values(["contig", "start"]).reset_index(drop=True)


def region_aln_df(
    db_path: str,
    assembly: str,
    contig: str,
    start: int,
    end: int,
    as_plot: bool = False,
    pid_threshold: float = 0.0,
    plen_threshold: float | None = None,
    pid_tie_eps: float = 0.0001,
) -> pd.DataFrame:
    """
    Convenience wrapper around segments_aln_df for a single contig region.
    Returns all alignments overlapping [start, end) on contig for assembly.

    as_plot : bool
        When True, apply the same post-processing as alignment_plot:
        1. Drop alignments below pid_threshold / plen_threshold.
        2. Per read, keep only alignments within pid_tie_eps of that read's
           best PID (ties get equal weight in the plot; default delta 0.0001).
    """
    from qatool.SegmentStorage import Segment, SegmentCollection
    sc = SegmentCollection([Segment(assembly, contig, start, end)])
    df = segments_aln_df(db_path, sc, assembly)

    if as_plot and not df.empty:
        df = df[df.pid >= pid_threshold]
        if plen_threshold is not None:
            df = df[df.plen >= plen_threshold]
        if not df.empty:
            best = df.groupby("read_name")["pid"].transform("max")
            df = df[df.pid + pid_tie_eps >= best]
        df = df.reset_index(drop=True)

    return df


def reads_all_aln_df(
    db_path: str,
    read_names,
    assembly: str | None = None,
) -> pd.DataFrame:
    """
    Return all alignments for the given reads across every assembly in the DB.

    Parameters
    ----------
    db_path : str
    read_names : iterable of str
    assembly : str or None
        If given, restrict to that assembly; otherwise returns all assemblies.

    Returns
    -------
    pd.DataFrame with columns: assembly, read_name, contig, start, end, strand,
        pid, plen, ed, q_aln_len, q_len, soft_left, soft_right, aln_type.
    Sorted by assembly, contig, start.
    """
    reads_list = list(read_names)
    if not reads_list:
        return pd.DataFrame(columns=["assembly"] + _ALN_COLS)

    rows_out = []
    con = sqlite3.connect(db_path)
    try:
        for i in range(0, len(reads_list), 800):
            chunk = reads_list[i:i + 800]
            ph = ",".join(["?"] * len(chunk))
            if assembly is not None:
                rows_out.extend(con.execute(
                    f"""SELECT assembly, read_name, contig, start, end, strand,
                               pid, plen, ed, q_aln_len, q_len,
                               soft_left, soft_right, aln_type
                        FROM aln
                        WHERE assembly = ? AND read_name IN ({ph})""",
                    [assembly] + chunk,
                ).fetchall())
            else:
                rows_out.extend(con.execute(
                    f"""SELECT assembly, read_name, contig, start, end, strand,
                               pid, plen, ed, q_aln_len, q_len,
                               soft_left, soft_right, aln_type
                        FROM aln
                        WHERE read_name IN ({ph})""",
                    chunk,
                ).fetchall())
    finally:
        con.close()

    df = pd.DataFrame(rows_out, columns=["assembly"] + _ALN_COLS)
    return df.sort_values(["assembly", "contig", "start"]).reset_index(drop=True)


def read_aln_df(
    db_path: str,
    read_name: str,
    assembly: str | None = None,
) -> pd.DataFrame:
    """
    Convenience wrapper around reads_all_aln_df for a single read.
    Returns all alignments for read_name, optionally filtered to one assembly.
    """
    return reads_all_aln_df(db_path, [read_name], assembly=assembly)


# ── sequence fetching ─────────────────────────────────────────────────────────

def fetch_seqs(
    fasta_path: str,
    regions,
) -> List[Tuple[str, int, str]]:
    """
    Fetch sequences from a bgz-indexed FASTA.

    Parameters
    ----------
    fasta_path : str
        Path to bgzip-compressed, samtools-indexed FASTA (.bgz).
    regions : iterable of (contig, start, end)

    Returns
    -------
    list of (contig, start, seq)
    """
    result = []
    fa = pysam.FastaFile(fasta_path)
    try:
        for contig, start, end in regions:
            result.append((str(contig), int(start), fa.fetch(str(contig), int(start), int(end))))
    finally:
        fa.close()
    return result

