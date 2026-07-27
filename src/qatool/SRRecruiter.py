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
import hashlib
import pysam
import pysam
import pandas as pd
import math
import os
from collections import Counter
from typing import Set, Iterable, List, Tuple, Dict
from qatool.SegmentStorage import *
from qatool.aln_utils import cigar_stats
from qatool.plotters import report_plot

class Recruiter:

    def __init__(self, db_path: str,
                 asm1_bam_path: str, asm2_bam_path: str,
                 asm1_fasta_path: str, asm2_fasta_path: str,
                 asm1: str = "asm1", asm2: str = "asm2",
                 merge_gap: int = 20000, verbose: bool = True,
                 pid_threshold: float=0.99, plen_threshold: float=0.99,
                 extension_window: int=100, extension_k: int=5,
                 ext_max: int=100_000, ext_rounds: int=100,
                 max_rounds=100, out_dir: str | None = "recruitment_out",
                 best_pid_delta: float = 0.0001):
        self.db_path = db_path
        self.con = sqlite3.connect(db_path)
        _tables = {r[0] for r in self.con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
        if "aln" not in _tables:
            raise RuntimeError(
                f"Database '{db_path}' has no 'aln' table. "
                "Run prepare_data.sh (or build_aln_db.py) to populate it first."
            )
        if asm1 == asm2:
            raise ValueError("asm1 and asm2 must be different")
        self.asm1 = asm1
        self.asm2 = asm2
        self.bams = {asm1: pysam.AlignmentFile(asm1_bam_path, "rb"),
                     asm2: pysam.AlignmentFile(asm2_bam_path, "rb")}
        self.fastas = {asm1: pysam.FastaFile(asm1_fasta_path),
                       asm2: pysam.FastaFile(asm2_fasta_path)}
        self.plen_threshold = plen_threshold
        self.pid_threshold = pid_threshold
        self.merge_gap = merge_gap
        self.max_rounds = max_rounds
        self.verbose = verbose
        self.extension_window = extension_window
        self.extension_k = extension_k
        self.ext_max = ext_max
        self.ext_rounds = ext_rounds
        self.out_dir = out_dir
        self.best_pid_delta = best_pid_delta

        self.sc: SegmentCollection = SegmentCollection(merge_gap=self.merge_gap)
        self.seen_reads: Set[str] = set()
        self.aln_df: pd.DataFrame = pd.DataFrame(columns=["read_name", "asm", "contig",
                                                          "start", "end", "pid", "plen", "ed"])
        self._sealed: Dict[Tuple[str, str, int, int], Tuple[bool, bool]] = {}
        self._contig_lens: Dict[str, Dict[str, int]] = {}
        self._common_cache: Dict[str, bool] = {}

    def reset(self):
        self.sc = SegmentCollection(merge_gap=self.merge_gap)
        self.seen_reads = set()
        self.aln_df = pd.DataFrame(columns=["read_name", "asm", "contig",
                                            "start", "end", "pid", "plen", "ed"])
        self._sealed = {}
        self._common_cache = {}

    def _concat_seqs_for_segments(self, sc_sorted: List[Segment], assembly: str) -> str:
        fa = self.fastas[assembly]
        parts: List[str] = []
        for seg in sc_sorted:
            if seg.assembly != assembly:
                continue
            parts.append(fa.fetch(seg.contig, seg.start, seg.end))
        return "".join(parts)

    def _write_round_plot(self, round_id: int):
        if self.out_dir is None:
            return
        os.makedirs(self.out_dir, exist_ok=True)
        sc_sorted = sorted(self.sc, key=lambda s: (s.contig, s.start, s.end))
        seq_asm1 = self._concat_seqs_for_segments(sc_sorted, self.asm1)
        seq_asm2 = self._concat_seqs_for_segments(sc_sorted, self.asm2)
        out_path = os.path.join(self.out_dir, f"round_{round_id}.png")
        """report_plot(
            seq_asm1, seq_asm2, self.sc, self.seen_reads,
            self.db_path, self.pid_threshold, self.plen_threshold,
            output_path=out_path, asm1=self.asm1, asm2=self.asm2, show=False
        )"""

    def run(self, seed: List[Segment]|SegmentCollection) -> Tuple[SegmentCollection, Set[str], pd.DataFrame]:
        """
        Runs segment and reads recruitment
        Output:
        segments - SegmentCollection object
        reads - Set of read names
        """
        self.reset()
        self.sc = SegmentCollection(seed, merge_gap=self.merge_gap)
        self._prime_contig_lengths()
        self._clamp_seed_to_contig_bounds()
        cur_round: int = 0
        while True:
            # assert(cur_round < self.max_rounds) # TODO! raise error if max rounds exceeded
            cur_round += 1
            if cur_round > self.max_rounds:
                print(f"max rounds ({self.max_rounds}) exseeded, stopping.")
                break
            if self.verbose: print(f"[round {cur_round} start]\tnumber of segments: {len(self.sc)}\t"+
                                                              f"number of seen reads: {len(self.seen_reads)}")
            if self.ext_rounds is None or cur_round <= self.ext_rounds:
                self.extension_phase()
            new_reads = self.expand_segments_to_reads()
            if not new_reads:
                if self.verbose: print(f"STOPPED no new reads")
                break
            self.seen_reads.update(new_reads)
            self.expand_reads_to_segments(new_reads)
            self._write_round_plot(cur_round)
        if self.verbose: print(f"DONE segments recruited:", self.sc)
        # Final pass: include reads fully contained in the final segments
        final_reads = self.expand_segments_to_reads()
        if final_reads:
            self.seen_reads.update(final_reads)
        return self.sc, self.seen_reads, self.aln_df

    def _prime_contig_lengths(self):
        if self._contig_lens:
            return
        for asm, bam in self.bams.items():
            self._contig_lens[asm] = {name: int(length) for name, length in zip(bam.references, bam.lengths)}

    def _clamp_seed_to_contig_bounds(self):
        clamped = []
        changed = False
        for seg in self.sc:
            contig_len = self._contig_lens.get(seg.assembly, {}).get(seg.contig)
            if contig_len is None:
                clamped.append(seg)
                continue
            new_start = max(0, seg.start)
            new_end = min(seg.end, contig_len)
            if new_start != seg.start or new_end != seg.end:
                print(f"WARNING: seed segment {seg.assembly}:{seg.contig}:{seg.start}-{seg.end} "
                      f"exceeds contig bounds [0, {contig_len}), clamped to {new_start}-{new_end}")
                changed = True
            clamped.append(Segment(seg.assembly, seg.contig, new_start, new_end))
        if changed:
            self.sc = SegmentCollection(clamped, merge_gap=self.merge_gap)

    def _remap_sealed_to_current_segments(self):
        if not self._sealed:
            return
        new_sealed: Dict[Tuple[str, str, int, int], Tuple[bool, bool]] = {}
        for seg in self.sc:
            key = (seg.assembly, seg.contig, seg.start, seg.end)
            if key in self._sealed:
                new_sealed[key] = self._sealed[key]
        self._sealed = new_sealed

    def _read_is_common(self, read_name: str) -> bool:
        cached = self._common_cache.get(read_name)
        if cached is not None:
            return cached

        rows = self.con.execute(
            """SELECT assembly, contig, start, end, pid
               FROM aln
               WHERE assembly IN (?,?)
                 AND read_name = ?
                 AND pid >= ?
                 AND plen >= ?""",
            [self.asm1, self.asm2, read_name, self.pid_threshold, self.plen_threshold],
        ).fetchall()

        # Step 1: best PID per assembly
        best_pid: Dict[str, float] = {}
        for asm, contig, start, end, pid in rows:
            pid = float(pid)
            if pid > best_pid.get(asm, float("-inf")):
                best_pid[asm] = pid

        if self.asm1 not in best_pid or self.asm2 not in best_pid:
            self._common_cache[read_name] = False
            return False

        if abs(best_pid[self.asm1] - best_pid[self.asm2]) > self.best_pid_delta:
            self._common_cache[read_name] = False
            return False

        # Step 2: collect best loci (within delta) per assembly
        loci: Dict[str, List[Tuple[str, int, int]]] = {self.asm1: [], self.asm2: []}
        for asm, contig, start, end, pid in rows:
            if float(pid) + self.best_pid_delta >= best_pid[asm]:
                loci[asm].append((str(contig), int(start), int(end)))

        # Step 3: compare sequence hash multisets
        _comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")

        def _seq_hashes(asm: str, loci_list: List[Tuple[str, int, int]]) -> Counter:
            fa = self.fastas[asm]
            out: Counter = Counter()
            for contig, start, end in loci_list:
                seq = fa.fetch(contig, start, end).upper()
                out[hashlib.sha1(seq.encode()).digest()] += 1
                out[hashlib.sha1(seq.translate(_comp)[::-1].encode()).digest()] += 1
            return out

        is_common = _seq_hashes(self.asm1, loci[self.asm1]) == _seq_hashes(self.asm2, loci[self.asm2])
        self._common_cache[read_name] = is_common
        return is_common

    def _check_window_for_extension(self, asm: str, contig: str, w_start: int, w_end: int) -> bool:
        if w_end <= w_start:
            return True
        common_reads: Set[str] = set()
        uncommon_reads: Set[str] = set()
        for aln in self.bams[asm].fetch(contig, w_start, w_end):
            assert(not aln.is_unmapped)
            soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid, _ = cigar_stats(aln)
            if pid < self.pid_threshold or plen is None or plen < self.plen_threshold:
                continue
            if self._read_is_common(aln.query_name):
                common_reads.add(aln.query_name)
            else:
                uncommon_reads.add(aln.query_name)
        if len(common_reads) > self.extension_k or not uncommon_reads:
            print("SEALED", asm, contig, w_start, w_end, len(common_reads))
            return True
        return False

    def _find_container_interval(self, sc: SegmentCollection, asm: str, contig: str,
                                 start: int, end: int) -> Tuple[int, int] | None:
        intervals = sc.storage.get((asm, contig), [])
        for s, e in intervals:
            if s <= start and e >= end:
                return s, e
        return None

    def extension_phase(self):
        self._remap_sealed_to_current_segments()
        if self.extension_window <= 0 or self.extension_k <= 0:
            return
        updated: List[Tuple[Segment, bool, bool]] = []

        for seg in list(self.sc):
            key = (seg.assembly, seg.contig, seg.start, seg.end)
            sealed_left, sealed_right = self._sealed.get(key, (False, False))
            start = seg.start
            end = seg.end
            orig_start = seg.start
            orig_end = seg.end
            contig_len = self._contig_lens[seg.assembly].get(seg.contig)
            if contig_len is None:
                raise ValueError(f"Unknown contig length for {seg.assembly}:{seg.contig}")
            if self.ext_max <= 0:
                left_limit = orig_start
                right_limit = orig_end
            else:
                left_limit = max(0, orig_start - self.ext_max)
                right_limit = min(contig_len, orig_end + self.ext_max)

            # extend left
            while not sealed_left:
                w_start = start
                w_end = min(start + self.extension_window, end)
                if start <= left_limit:
                    sealed_left = True
                    break
                if w_end <= w_start:
                    sealed_left = True
                    break
                left_ends: List[int] = []
                common_left_ends: List[int] = []
                for aln in self.bams[seg.assembly].fetch(seg.contig, w_start, w_end):
                    _, _, _, _, _, _, plen, pid, _ = cigar_stats(aln)
                    if pid < self.pid_threshold or plen is None or plen < self.plen_threshold:
                        continue
                    left_ends.append(int(aln.reference_start))
                    if self._read_is_common(aln.query_name):
                        common_left_ends.append(int(aln.reference_start))
                if not left_ends:
                    sealed_left = True
                    break
                if self._check_window_for_extension(seg.assembly, seg.contig, w_start, w_end):
                    # Window is common-dominated: snap out to the nearest bridging common
                    # read (the one whose start is closest to the boundary), then seal.
                    # Do not keep crawling through the common region.
                    candidates = [s for s in common_left_ends if s <= start]
                    if candidates:
                        new_start = max(max(candidates), left_limit)
                        if new_start != start:
                            print("LEFT_ext", seg.assembly, seg.contig, w_start, w_end)
                            start = new_start
                    sealed_left = True
                    break
                # Discriminating window: extend to the outermost qualifying read.
                new_start = max(left_limit, min(left_ends))
                if new_start == start:
                    sealed_left = True
                    break
                print("LEFT_ext", seg.assembly, seg.contig, w_start, w_end)
                start = new_start

            # extend right
            while not sealed_right:
                w_end = end
                w_start = max(end - self.extension_window, start)
                if end >= right_limit:
                    sealed_right = True
                    break
                if w_end <= w_start:
                    sealed_right = True
                    break
                right_ends: List[int] = []
                common_right_ends: List[int] = []
                for aln in self.bams[seg.assembly].fetch(seg.contig, w_start, w_end):
                    _, _, _, _, _, _, plen, pid, _ = cigar_stats(aln)
                    if pid < self.pid_threshold or plen is None or plen < self.plen_threshold:
                        continue
                    right_ends.append(int(aln.reference_end))
                    if self._read_is_common(aln.query_name):
                        common_right_ends.append(int(aln.reference_end))
                if not right_ends:
                    sealed_right = True
                    break
                if self._check_window_for_extension(seg.assembly, seg.contig, w_start, w_end):
                    # Window is common-dominated: snap out to the nearest bridging common
                    # read (the one whose end is closest to the boundary), then seal.
                    # Do not keep crawling through the common region.
                    candidates = [e for e in common_right_ends if e >= end]
                    if candidates:
                        new_end = min(min(candidates), right_limit)
                        if new_end != end:
                            print("RIGHT_ext", seg.assembly, seg.contig, w_start, w_end)
                            end = new_end
                    sealed_right = True
                    break
                # Discriminating window: extend to the outermost qualifying read.
                new_end = min(right_limit, max(right_ends))
                if new_end == end:
                    sealed_right = True
                    break
                print("RIGHT_ext", seg.assembly, seg.contig, w_start, w_end)
                end = new_end

            updated.append((Segment(seg.assembly, seg.contig, start, end), sealed_left, sealed_right))

        new_sc = SegmentCollection(merge_gap=self.merge_gap)
        for seg, _, _ in updated:
            new_sc.add(seg)
        new_sealed: Dict[Tuple[str, str, int, int], Tuple[bool, bool]] = {}
        for seg in new_sc:
            new_sealed[(seg.assembly, seg.contig, seg.start, seg.end)] = (False, False)

        for seg, sealed_left, sealed_right in updated:
            container = self._find_container_interval(new_sc, seg.assembly, seg.contig, seg.start, seg.end)
            if container is None:
                continue
            c_start, c_end = container
            key = (seg.assembly, seg.contig, c_start, c_end)
            cur_left, cur_right = new_sealed.get(key, (False, False))
            if sealed_left and seg.start == c_start:
                cur_left = True
            if sealed_right and seg.end == c_end:
                cur_right = True
            new_sealed[key] = (cur_left, cur_right)

        self.sc = new_sc
        self._sealed = new_sealed

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
                if aln.query_name in self.seen_reads:
                    continue
                #pid = calculate_pid(aln)
                soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid, _ = cigar_stats(aln)
                if pid < self.pid_threshold or plen is None or plen < self.plen_threshold:
                    continue
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
        raise ValueError("assembly_lengths must be provided, e.g. {'asm1': L_asm1, 'asm2': L_asm2}")

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
