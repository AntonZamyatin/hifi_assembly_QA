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
from typing import Set, Iterable
from SegmentStorage import *
from aln_utils import calculate_pid, cigar_stats

class Recruiter:

    def __init__(self,  seed: List[Segment], db_path: str,
                 ref_bam_path: str, cand_bam_path: str,
                 merge_gap: int, verbose: bool):
        self.con = sqlite3.connect(db_path)
        self.seed = seed
        self.bams = {"ref": pysam.AlignmentFile(ref_bam_path, "rb"),
                     "cand": pysam.AlignmentFile(cand_bam_path, "rb")}
        self.plen_threshold = 0.9
        self.pid_primary = 0.99
        self.pid_secondary = 0.92
        self.pid_supp = 1
        self.min_supp_q_aln_len = 1000
        self.merge_gap = merge_gap
        self.max_rounds = 100
        self.verbose = verbose

        self.sc: SegmentCollection = SegmentCollection(merge_gap=self.merge_gap)
        self.seen_reads: Set[str] = set()
        self.aln_df: pd.DataFrame = pd.DataFrame(columns=["read_name", "asm", "contig",
                                                          "start", "end", "pid", "plen"])

    def reset(self):
        self.sc = SegmentCollection(merge_gap=self.merge_gap)
        self.seen_reads = set()
        self.aln_df = pd.DataFrame(columns=["read_name", "asm", "contig",
                                            "start", "end", "pid", "plen"])

    def run(self) -> Tuple[SegmentCollection, Set[str], pd.DataFrame]:
        """
        Runs segment and reads recruitment
        Output:
        segments - SegmentCollection object
        reads - Set of read names
        """
        self.reset()
        self.sc = SegmentCollection(self.seed, merge_gap=self.merge_gap)
        cur_round: int = 0
        while cur_round != self.max_rounds:
            cur_round += 1
            if self.verbose: print(f"[round {cur_round} start]\tnumber of segments: {len(self.sc)}\t"+
                                                              f"number of seen reads: {len(self.seen_reads)}")
            new_reads = self.expand_segments_to_reads()
            if not new_reads:
                if self.verbose: print(f"STOPPED no new reads")
                break
            self.seen_reads.update(new_reads)
            self.expand_reads_to_segments(new_reads)
            if self.verbose: print(f"[round {cur_round} start]\tnumber of segments: {len(self.sc)}\t"+
                                                              f"number of seen reads: {len(self.seen_reads)}")
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
                if aln.is_unmapped or aln.reference_start < seg.start or aln.reference_end > seg.end:
                    continue
                # secondary alignments:
                pid = calculate_pid(aln)
                soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen = cigar_stats(aln)
                if aln.is_secondary:
                    if pid < self.pid_secondary or plen is None or plen < self.plen_threshold:
                        continue
                    if aln.query_name not in self.seen_reads:
                        reads.add(aln.query_name)
                        continue
                # supplementary alignments
                if aln.is_supplementary:
                    if pid < self.pid_supp or q_aln_len < self.min_supp_q_aln_len:
                        continue
                    if aln.query_name not in self.seen_reads:
                        reads.add(aln.query_name)
                        continue
                # primary alignments
                if pid < self.pid_primary or plen is None or plen < self.plen_threshold:
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
                    SELECT read_name, assembly, contig, start, end, pid, plen
                    FROM aln
                    WHERE assembly IN ({asm_placeholders})
                      AND contig IS NOT NULL
                      AND start IS NOT NULL
                      AND end IS NOT NULL
                      AND read_name IN ({{READS_IN}})
                      AND (
                            (aln_type IN ('primary','primary_wsup') AND pid >= ? AND plen >= ?)
                         OR (aln_type = 'secondary'               AND pid >= ? AND plen >= ?)
                         OR (aln_type = 'supplementary'           AND pid >= ? AND q_aln_len >= ?)
                      )
                    """
        tail_params = [self.pid_primary, self.plen_threshold,
                       self.pid_secondary, self.plen_threshold,
                       self.pid_supp, self.min_supp_q_aln_len]
        for chunk in _chunks(list(reads), 1000):
            reads_in = ",".join(["?"] * len(chunk))
            sql = base_sql.replace("{READS_IN}", reads_in)
            params = list(self.bams.keys()) + chunk + tail_params
            rows = self.con.execute(sql, params).fetchall()
            for read_name, asm, contig, start, end, pid, plen in rows:
                new_segments.add(Segment(str(asm), str(contig), int(start), int(end)))
                aln_rows.append({"read_name": str(read_name), "asm": str(asm),
                                 "contig": str(contig), "start": int(start),
                                 "end": int(end), "pid": float(pid), "plen": float(plen)})
        self.aln_df = pd.concat([self.aln_df, pd.DataFrame(aln_rows)], ignore_index=True)
        self.sc += new_segments
