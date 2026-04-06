"""
Project: hifi assembly QA
Module: test.py
Description: Implements recursive segment–read recruitment between reference
             and candidate assemblies using HiFi read alignments stored in
             indexed BAM files and an accompanying SQLite alignment database.
             The module iteratively expands genomic segments via supporting
             reads and projects them across assemblies until convergence.

Author: Anton Zamyatin
Affiliation: Pennsylvania State University
Contact: zamyatin.anton@gmail.com

Created: 2/11/26

License: MIT License
"""

import sqlite3
import pysam
import pandas as pd
import math
from typing import Set, Iterable, List, Tuple, Dict
from SegmentStorage import *
from aln_utils import cigar_stats

class Recruiter:

    def __init__(self,  db_path: str,
                 ref_bam_path: str, cand_bam_path: str,
                 merge_gap: int, verbose: bool,
                 pid_threshold: float=0.99, plen_threshold: float=0.99):
        self.con = sqlite3.connect(db_path)
        self.bams = {"ref": pysam.AlignmentFile(ref_bam_path, "rb"),
                     "cand": pysam.AlignmentFile(cand_bam_path, "rb")}
        self.plen_threshold = plen_threshold
        self.pid_threshold = pid_threshold
        self.merge_gap = merge_gap
        self.max_rounds = 100
        self.verbose = verbose

        self.sc: SegmentCollection = SegmentCollection(merge_gap=self.merge_gap)
        self.seen_reads: Set[str] = set()
        self.aln_df: pd.DataFrame = pd.DataFrame(columns=["read_name", "asm", "contig",
                                                          "start", "end", "pid", "plen", "ed"])

    def reset(self):
        self.sc = SegmentCollection(merge_gap=self.merge_gap)
        self.seen_reads = set()
        self.aln_df = pd.DataFrame(columns=["read_name", "asm", "contig",
                                            "start", "end", "pid", "plen", "ed"])

    def run(self, seed: List[Segment]|SegmentCollection) -> Tuple[SegmentCollection, Set[str], pd.DataFrame]:
        """
        Runs segment and reads recruitment
        Output:
        segments - SegmentCollection object
        reads - Set of read names
        """
        self.reset()
        self.sc = SegmentCollection(seed, merge_gap=self.merge_gap)
        cur_round: int = 0
        while True:
            assert(cur_round < self.max_rounds) # TODO! raise error if max rounds exceeded
            cur_round += 1
            if self.verbose: print(f"[round {cur_round} start]\tnumber of segments: {len(self.sc)}\t"+
                                                              f"number of seen reads: {len(self.seen_reads)}")
            new_reads = self.expand_segments_to_reads()
            if not new_reads:
                if self.verbose: print(f"STOPPED no new reads")
                break
            self.seen_reads.update(new_reads)
            self.expand_reads_to_segments(new_reads)
        if self.verbose: print(f"DONE segments recruited:", self.sc)
        return self.sc, self.seen_reads, self.aln_df

    def expand_segments_to_reads(self) -> Set[str]:
        """
        Fetch read names from the alignments that are contained in segments boundaries and satisfy
        identity and length criteria.
        Return:
            Set of read names
        """
        reads: Set[str] = set()
        for seg in self.sc:
            for aln in self.bams[seg.assembly].fetch(seg.contig, seg.start, seg.end):
                assert(not aln.is_unmapped)
                if aln.reference_start < seg.start or aln.reference_end > seg.end:
                    continue
                #pid = calculate_pid(aln)
                soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid = cigar_stats(aln)
                if pid < self.pid_threshold or plen is None or plen < self.plen_threshold:
                    continue
                if aln.query_name not in self.seen_reads:
                    reads.add(aln.query_name)
        return reads

    def expand_reads_to_segments(self, reads: Set[str]):
        """
        Fetch segments covered by reads from the alignment database indexed by read names.
        Returns:
            SegmentCollection
        """
        def _chunks(lst: List[str], n: int) -> Iterable[List[str]]:
            for i in range(0, len(lst), n):
                yield lst[i:i + n]

        new_segments: SegmentCollection = SegmentCollection()
        aln_rows = []
        asm_placeholders: str = ",".join(["?"] * len(self.bams.keys()))
        base_sql = f"""
                    SELECT read_name, assembly, contig, start, end, pid, plen, ed
                    FROM aln
                    WHERE assembly IN ({asm_placeholders})
                      AND contig IS NOT NULL
                      AND start IS NOT NULL
                      AND end IS NOT NULL
                      AND read_name IN ({{READS_IN}})
                      AND pid >= ?
                      AND plen >= ?
                    """
        tail_params = [self.pid_threshold, self.plen_threshold]
        for chunk in _chunks(list(reads), 1000):
            reads_in = ",".join(["?"] * len(chunk))
            sql = base_sql.replace("{READS_IN}", reads_in)
            params = list(self.bams.keys()) + chunk + tail_params
            rows = self.con.execute(sql, params).fetchall()
            for read_name, asm, contig, start, end, pid, plen, ed in rows:
                new_segments.add(Segment(str(asm), str(contig), int(start), int(end)))
                aln_rows.append({"read_name": str(read_name), "asm": str(asm),
                                 "contig": str(contig), "start": int(start),
                                 "end": int(end), "pid": float(pid), "plen": float(plen)})
                aln_rows[-1]["ed"] = (None if ed is None else int(ed))
        self.aln_df = pd.concat([self.aln_df, pd.DataFrame(aln_rows)], ignore_index=True)
        self.sc += new_segments


def metrics_report_from_aln_df(
    aln_df: pd.DataFrame,
    assembly_lengths: Dict[str, float],
    pid_threshold: float | None = None,
    plen_threshold: float | None = None,
) -> pd.DataFrame:
    """
    Compute per-assembly support/likelihood metrics from recruiter alignment table.

    Metrics per assembly:
      - reads_only: number of reads with qualifying alignments only in that assembly
      - likelihood: Π_r ( (Σ_i exp(-D_i)) / L_A )
      - log_likelihood: Σ_r log( (Σ_i exp(-D_i)) / L_A )

    Uses D_i = ed if 'ed' exists; otherwise falls back to D_i = (1 - pid) * q_aln_len
    when 'q_aln_len' is present.
    """
    required = {"read_name", "asm"}
    missing = required - set(aln_df.columns)
    if missing:
        raise ValueError(f"aln_df is missing required columns: {sorted(missing)}")
    if not assembly_lengths:
        raise ValueError("assembly_lengths must be provided, e.g. {'ref': L_ref, 'cand': L_cand}")

    df = aln_df.copy()
    if pid_threshold is not None and "pid" in df.columns:
        df = df[df["pid"] >= float(pid_threshold)]
    if plen_threshold is not None and "plen" in df.columns:
        df = df[df["plen"] >= float(plen_threshold)]

    if "ed" in df.columns:
        df = df[df["ed"].notna()].copy()
        df["D"] = df["ed"].astype(float)
    elif "q_aln_len" in df.columns and "pid" in df.columns:
        print("here")
        df["D"] = (1.0 - df["pid"].astype(float)) * df["q_aln_len"].astype(float)
    else:
        raise ValueError("Need either 'ed' column or both 'pid' and 'q_aln_len' to compute D")

    reads = sorted(df["read_name"].astype(str).unique().tolist())
    if not reads:
        return pd.DataFrame(columns=["asm", "reads_only", "likelihood", "log_likelihood", "n_reads", "L_A"])

    read_sets: Dict[str, set] = {}
    for asm, g in df.groupby("asm"):
        read_sets[str(asm)] = set(g["read_name"].astype(str).unique().tolist())

    rows = []
    all_assemblies = sorted(read_sets.keys())
    for asm in all_assemblies:
        l_a = float(assembly_lengths.get(asm, 0.0))
        if l_a <= 0:
            raise ValueError(f"assembly_lengths['{asm}'] must be > 0")

        own_reads = read_sets.get(asm, set())
        other_reads = set()
        for other_asm, rr in read_sets.items():
            if other_asm != asm:
                other_reads |= rr
        reads_only = len(own_reads - other_reads)

        by_read_sum = (
            df[df["asm"].astype(str) == asm]
            .assign(_exp=lambda x: x["pid"].map(lambda d: math.exp(-float(d))))
            .groupby("read_name", as_index=True)["_exp"]
            .sum()
            .to_dict()
        )

        log_likelihood = 0.0
        valid = True
        for r in reads:
            p = float(by_read_sum.get(r, 0.0)) / l_a
            if p <= 0.0:
                continue
            log_likelihood += math.log(p)
        likelihood = 0.0 if not valid else math.exp(log_likelihood)
        if not valid:
            log_likelihood = float("-inf")

        rows.append({
            "asm": asm,
            "reads_only": int(reads_only),
            "likelihood": float(likelihood),
            "log_likelihood": float(log_likelihood),
            "n_reads": int(len(reads)),
            "L_A": float(l_a),
        })

    return pd.DataFrame(rows).sort_values("asm").reset_index(drop=True)
