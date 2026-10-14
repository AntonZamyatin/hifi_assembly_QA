"""Synthetic BAM / SQLite fixtures for the compute-layer tests.

Everything is generated in a temp dir — the tests never touch the real alignment
DBs, BAMs or FASTAs.
"""

import math
import os
import sqlite3

import pysam

from qatool.SegmentStorage import Segment

CONTIGS = {"ctgA": 20_000, "ctgB": 12_000}

BINSIZE = 1000
BEST_PID_DELTA = 0.0001

# Segments deliberately not multiples of BINSIZE, so the last bin of each is partial.
SEGS_A = [Segment("asm1", "ctgA", 2_000, 6_500), Segment("asm1", "ctgA", 9_000, 11_200)]
SEGS_B = [Segment("asm2", "ctgB", 1_000, 4_300)]

# (name, contig, start, length, nm, flag);  flag 256 = secondary, 2048 = supplementary
READS_A = [
    ("r1", "ctgA", 2_000, 1_500, 0, 0),
    ("r2", "ctgA", 2_500, 3_000, 2, 0),
    ("r3", "ctgA", 5_000, 1_400, 0, 0),
    ("r4", "ctgA", 9_100, 2_000, 1, 0),
    ("r5", "ctgA", 3_000, 1_000, 0, 0),
    ("r5", "ctgA", 9_500, 1_000, 0, 256),      # multi-mapper: weight 0.5 each
    ("r6", "ctgA", 4_000, 800, 400, 0),        # pid 0.5 -> filtered out
]
READS_B = [
    ("r1", "ctgB", 1_000, 1_500, 0, 0),
    ("r2", "ctgB", 1_800, 2_000, 0, 0),
    ("r3", "ctgB", 3_000, 1_200, 0, 0),
    ("r7", "ctgB", 1_200, 900, 0, 0),
]


def write_bam(path: str, reads: list, contigs: dict = CONTIGS) -> str:
    """Write and index a BAM.

    reads: list of (name, contig, start, length, nm, flag[, extra]).  CIGAR is a
    plain `<length>M`, so cigar_stats takes the NM branch and pid = 1 - nm/length.
    `extra` (optional dict): "clip": (left, right) hard-clipped bases, for a piece of
    a split read (plen = length / (length + clips)); "tags": more tags, e.g. SA.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": ln} for c, ln in contigs.items()]}
    order = {c: i for i, c in enumerate(contigs)}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for name, contig, start, length, nm, flag, *rest in sorted(
                reads, key=lambda r: (order[r[1]], r[2])):
            extra = rest[0] if rest else {}
            left, right = extra.get("clip", (0, 0))
            a = pysam.AlignedSegment()
            a.query_name = name
            a.reference_id = order[contig]
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = ([(5, left)] if left else []) + [(0, length)] + ([(5, right)] if right else [])
            a.flag = flag
            a.query_sequence = "A" * length
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            a.set_tag("NM", nm)
            for k, v in extra.get("tags", {}).items():
                a.set_tag(k, v)
            out.write(a)
    pysam.index(path)
    return path


def write_db(path: str, rows: list) -> str:
    """Build a minimal `aln` table.

    rows: list of dicts with at least assembly, read_name, contig, start, end.
    Remaining columns take sensible defaults (perfect, forward, primary).
    """
    con = sqlite3.connect(path)
    con.execute("""
        CREATE TABLE aln (
            assembly TEXT NOT NULL, read_name TEXT NOT NULL, contig TEXT,
            start INTEGER, end INTEGER, strand INTEGER, cigar TEXT, pid REAL,
            soft_left INTEGER NOT NULL, soft_right INTEGER NOT NULL,
            hard_left INTEGER NOT NULL, hard_right INTEGER NOT NULL,
            q_aln_len INTEGER NOT NULL, q_len INTEGER NOT NULL,
            ed INTEGER, plen REAL, aln_type TEXT NOT NULL
        )""")
    for r in rows:
        span = r["end"] - r["start"]
        d = {"strand": 1, "cigar": f"{span}M", "pid": 1.0,
             "soft_left": 0, "soft_right": 0, "hard_left": 0, "hard_right": 0,
             "q_aln_len": span, "q_len": span, "ed": 0, "plen": 1.0,
             "aln_type": "primary", **r}
        con.execute(
            "INSERT INTO aln VALUES (:assembly,:read_name,:contig,:start,:end,:strand,"
            ":cigar,:pid,:soft_left,:soft_right,:hard_left,:hard_right,:q_aln_len,"
            ":q_len,:ed,:plen,:aln_type)", d)
    con.commit()
    con.close()
    return path


def make_region_data(tmp: str, prefix: str = "rd", binsize: int = BINSIZE):
    """A small but structurally complete RegionData built entirely from fixtures."""
    from qaview.region import RegionData

    bam1 = write_bam(os.path.join(tmp, f"{prefix}1.bam"), READS_A)
    bam2 = write_bam(os.path.join(tmp, f"{prefix}2.bam"), READS_B)
    rows = []
    for name, contig, start, length, _nm, flag in READS_A:
        rows.append({"assembly": "asm1", "read_name": name, "contig": contig,
                     "start": start, "end": start + length,
                     "aln_type": "secondary" if flag & 256 else "primary"})
    for name, contig, start, length, _nm, _flag in READS_B:
        rows.append({"assembly": "asm2", "read_name": name, "contig": contig,
                     "start": start, "end": start + length, "aln_type": "primary"})
    db = write_db(os.path.join(tmp, f"{prefix}.sqlite"), rows)

    return RegionData.build(
        SEGS_A + SEGS_B, {"r1", "r2", "r3"}, db, bam1, bam2,
        asm1="asm1", asm2="asm2", pid_threshold=0.99, plen_threshold=0.5,
        binsize=binsize, best_pid_delta=BEST_PID_DELTA)


def marginal(rd, asm: str, bands=None, n_bands: int = 1, binsize: int | None = None) -> dict:
    """`asm`'s marginal over its whole axis at `binsize` (default rd.binsize), as the
    page's first frame bins it (region.marginal_bins); one band unless given."""
    from qaview.region import marginal_bins
    ext = rd.plot_h if asm == rd.asm1 else rd.plot_w
    bands = [0] * len(rd.read_names) if bands is None else bands
    return marginal_bins(rd.segs_info(asm), rd.placements[asm], bands, n_bands, 0, ext,
                         binsize or rd.binsize)


def reference_coverage_bins(bamfile, segments, pid_threshold, plen_threshold,
                            binsize, read_buckets, best_pid_delta):
    """Pre-refactor `alignment_plot` binning, verbatim, as a regression oracle.

    Copied from git 6c9e776 `plotters/coverage.py`.  Returns (bins, bucket_bins).
    """
    from qatool.aln_utils import cigar_stats

    segments = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    bins = [0.0 for _ in range(sum(math.ceil((s.end - s.start) / binsize) for s in segments))]
    bucket_bins = None
    name_to_bucket = None
    if read_buckets is not None:
        bucket_bins = [[0.0 for _ in range(len(bins))] for _ in range(len(read_buckets))]
        name_to_bucket = {}
        for i, bucket in enumerate(read_buckets):
            for name in bucket:
                name_to_bucket[name] = i

    valid_alns = []
    offset = 0
    bam = pysam.AlignmentFile(bamfile, "rb")
    try:
        for seg in segments:
            for aln in bam.fetch(seg.contig, start=seg.start, end=seg.end):
                _, _, _, _, q_aln_len, q_len, plen, pid, _ = cigar_stats(aln)
                if pid < pid_threshold:
                    continue
                if plen_threshold is not None and (plen is None or plen < plen_threshold):
                    continue
                startbin = math.floor((max(aln.reference_start, seg.start) - seg.start + offset) / binsize)
                endbin = math.ceil((min(aln.reference_end, seg.end) - seg.start + offset) / binsize)
                if startbin < endbin:
                    valid_alns.append((aln.query_name, float(pid), startbin, endbin))
            offset += math.ceil((seg.end - seg.start) / binsize) * binsize
    finally:
        bam.close()

    best_pid_by_read = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if read_name not in best_pid_by_read or pid > best_pid_by_read[read_name]:
            best_pid_by_read[read_name] = pid

    qualifying = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if pid + best_pid_delta >= best_pid_by_read.get(read_name, float("-inf")):
            qualifying.setdefault(read_name, []).append((startbin, endbin))

    for read_name, ranges in qualifying.items():
        w = 1.0 / len(ranges)
        if bucket_bins is None:
            for startbin, endbin in ranges:
                for i in range(startbin, endbin):
                    bins[i] += w
        else:
            bucket_idx = name_to_bucket.get(read_name)
            target = bucket_bins[bucket_idx] if bucket_idx is not None else bins
            for startbin, endbin in ranges:
                for i in range(startbin, endbin):
                    target[i] += w

    return bins, bucket_bins
