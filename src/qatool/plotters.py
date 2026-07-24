"""
Project: hifi assembly QA
Module: plotters.py
Description:

Author: Anton Zamyatin
Affiliation: Pennsylvania State University
Contact: zamyatin.anton@gmail.com

Created: 2/25/26

License: MIT License
"""

import hashlib
import pysam
import matplotlib.pyplot as plt
import math
import numpy as np
import sqlite3
import matplotlib.ticker as mticker
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory, offset_copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from qatool.aln_utils import cigar_stats
from typing import List
from qatool.SegmentStorage import SegmentCollection, Segment


PID_TIE_EPS_DEFAULT = 0.001
BEST_PID_DELTA_DEFAULT = 0.0001

SUPPORT_COLORS = {
    "asm1_only": "#b2182b",   # strong red
    "asm1_best": "#f4a3c2",   # light pink
    "asm2_only": "#1a9641",   # strong green
    "asm2_best": "#a6d96a",   # light green
    "common":    "#a3a395",   # grey-green
    "other":     "#424241",   # dark grey
}

# Occurrence classification colors: indices 0-4 = 1/2/3/4/5+ best loci (non-common),
# indices 5-9 = same counts for common reads (pale versions).
OCCURRENCE_COLORS = [
    "#2166ac",  # 1 locus
    "#92c5de",  # 1 locus, common
    "#1a9641",  # 2 loci
    "#a6d96a",  # 2 loci, common
    "#d7191c",  # 3 loci
    "#f4a582",  # 3 loci, common
    "#e6ab02",  # 4 loci
    "#fee08b",  # 4 loci, common
    "#252525",  # 5+ loci
    "#969696",  # 5+ loci, common
]


def _build_axis_map(sc, assembly: str, binsize: int = None):
    """
    Build a (contig, start, end, offset) list for segments of the given assembly,
    in the order they appear in sc. Returns the list and the total length.

    If binsize is given, each segment's contribution is rounded UP to the next
    bin boundary (ceil(seg_len / binsize) * binsize), so all offsets are exact
    multiples of binsize and segment boundaries coincide with bin edges.
    """
    segs = []
    offset = 0
    for seg in sc:
        if seg.assembly != assembly:
            continue
        segs.append((seg.contig, seg.start, seg.end, offset))
        seg_len = seg.end - seg.start
        if binsize is not None:
            offset += math.ceil(seg_len / binsize) * binsize
        else:
            offset += seg_len
    return segs, offset


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


def alignment_plot(bamfile: str, segments: List[Segment],
                   pid_threshold: float = 0.99, plen_threshold: float | None = None,
                   binsize: int = 1000,
                   ax=None, read_buckets=None, read_bucket_colors=None, orientation: str = "x",
                   return_max=False, best_pid_delta: float = BEST_PID_DELTA_DEFAULT):
    """
    Generate a simple alignment plot for a given BAM file, contig, and region.

    Parameters:
    - bamfile: Path to the BAM file containing alignments.
    - segments: List of Segment to plot alignments for (same assembly/contig).
    """
    if not segments:
        raise ValueError("segments must be a non-empty list")
    segments = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    # Each segment occupies ceil(seg_len/binsize) bins so that segment offsets
    # land exactly on bin boundaries (no bin straddles a segment edge).
    bins = [0.0 for _ in range(sum(math.ceil((s.end - s.start) / binsize) for s in segments))]
    bucket_bins = None
    name_to_bucket = None
    if read_buckets is not None:
        bucket_bins = [[0.0 for _ in range(len(bins))] for _ in range(len(read_buckets))]
        name_to_bucket = {}
        for i, bucket in enumerate(read_buckets):
            for name in bucket:
                name_to_bucket[name] = i

    # Collect valid alignments; later keep only best-PID per read.
    valid_alns = []
    offset = 0
    bam = pysam.AlignmentFile(bamfile, "rb")
    try:
        for seg in segments:
            for aln in bam.fetch(seg.contig, start=seg.start, end=seg.end):
                soft_left, soft_right, hard_left, hard_right, q_aln_len, q_len, plen, pid, _ = cigar_stats(aln)
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

    # Pass 1: find best pid per read.
    best_pid_by_read = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if read_name not in best_pid_by_read or pid > best_pid_by_read[read_name]:
            best_pid_by_read[read_name] = pid

    # Pass 2: collect qualifying alignments (within pid_tie_eps of best) per read,
    # then fill bins with weight 1 / number of tied alignments.
    qualifying: dict = {}
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

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(16, 5))

    # n+1 x-values so step='post' renders all n bins fully:
    # the last bin covers [x[n-1], x[n]] and needs a right edge.
    x = np.arange(len(bins) + 1) * binsize
    if bucket_bins is None:
        bins_plot = bins + [0.0]
        if orientation == "x":
            ax.stackplot(x, bins_plot, color='blue', step='post')
        elif orientation == "y":
            ax.fill_betweenx(x, 0, bins_plot, step='post', color='blue', linewidth=0)
        else:
            raise ValueError("orientation must be 'x' or 'y'")
    else:
        if read_bucket_colors is None:
            colors = [f"C{i % 10}" for i in range(len(bucket_bins))] + ['#9e9e9e']
        else:
            if len(read_bucket_colors) != len(bucket_bins):
                raise ValueError("read_bucket_colors length must match read_buckets length")
            colors = list(read_bucket_colors) + ['#6e6e6e']
        series = [s + [0.0] for s in bucket_bins + [bins]]
        if orientation == "x":
            ax.stackplot(x, series, colors=colors, step='post')
        elif orientation == "y":
            y0 = np.zeros(len(x))
            for s, c in zip(series, colors):
                y1 = y0 + np.array(s)
                ax.fill_betweenx(x, y0, y1, step='post', color=c, linewidth=0)
                y0 = y1
        else:
            raise ValueError("orientation must be 'x' or 'y'")

    # Axes styling: no rectangle, axes through origin
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    if standalone:
        if orientation == "x":
            ax.set_xlim(0, max(len(bins) * binsize - 1, 0))
        else:
            ax.set_ylim(0, max(len(bins) * binsize - 1, 0))
        plt.xlabel('Position in Bin')
        plt.ylabel('Number of Alignments')
        plt.title('Alignment Plot for Segments (' + str(len(segments)) + ')')
        plt.show()

    if return_max:
        if bucket_bins is None:
            return ax, max(bins) if bins else 0.0, None
        total = np.array(bins, dtype=float)
        for s in bucket_bins:
            total += np.array(s, dtype=float)
        return ax, float(np.max(total)) if len(total) > 0 else 0.0, bucket_bins
    return ax


def dotplot_kmer(seq_a: str, seq_b: str, k: int = 15, binsize: int = 1000,
                 ax=None, max_cells: int = 10_000_000,
                 mask_ambiguous: bool = True, log_scale: bool = False,
                 max_kmer_freq: int = 100):
    """
    Generate a Gepard-like k-mer dot plot using binned k-mer matches.

    Parameters:
    - seq_a: First sequence string.
    - seq_b: Second sequence string.
    - k: k-mer size.
    - binsize: Window size for binning positions.
    - max_cells: Max allowed dense matrix cells (bins_a * bins_b).
    - max_kmer_freq: k-mers appearing more than this many times in seq_a are
        excluded from the index (repeat / low-complexity filter). 0 = no filter.
        Default 100 is appropriate for genome-scale sequences; increase if the
        plot is missing signal in repetitive regions.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    if binsize <= 0:
        raise ValueError("binsize must be positive")

    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    len_a = len(seq_a)
    len_b = len(seq_b)
    if len_a < k or len_b < k:
        raise ValueError("Sequences must be at least k long")

    bins_a = math.ceil(len_a / binsize)
    bins_b = math.ceil(len_b / binsize)
    if bins_a * bins_b > max_cells:
        raise ValueError(
            f"dotplot matrix would be {bins_a}x{bins_b} cells; "
            "increase binsize or max_cells"
        )

    # Byte-level lookup: ASCII -> 2-bit base (A=0,C=1,G=2,T=3), ambiguous -> 255.
    enc_table = np.full(256, 255, dtype=np.uint8)
    for _b, _v in zip(b'ACGT', [0, 1, 2, 3]):
        enc_table[_b] = _v

    def _kmer_hashes(enc: np.ndarray, do_rc: bool) -> tuple:
        """
        Compute k-mer hashes for all valid positions using a k-iteration vectorised loop.
        Each iteration adds one shifted base, consuming only O(n) memory (no (n,k) array).
        Returns (positions, fwd_hashes, rc_hashes|None).
        """
        n = len(enc)
        n_pos = n - k + 1
        if n_pos <= 0:
            empty = np.array([], dtype=np.int64)
            return np.array([], dtype=np.intp), empty, empty if do_rc else None

        ph = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)  # fwd: 4^(k-1)..4^0
        pl = 4 ** np.arange(k,           dtype=np.int64)     # rc:  4^0..4^(k-1)
        ei = enc.astype(np.int64)

        fwd = np.zeros(n_pos, dtype=np.int64)
        rci = np.zeros(n_pos, dtype=np.int64) if do_rc else None
        for j in range(k):
            sl = ei[j: j + n_pos]
            fwd += sl * ph[j]
            if do_rc:
                rci += sl * pl[j]

        if mask_ambiguous:
            # Mark windows that contain any ambiguous base (255 in encoded form).
            inv_pos = np.where(enc == 255)[0]
            valid = np.ones(n_pos, dtype=bool)
            for p in inv_pos:
                valid[max(0, p - k + 1): min(n_pos, p + 1)] = False
            pos = np.where(valid)[0]
            fwd = fwd[valid]
            rc  = (np.int64(4 ** k - 1) - rci[valid]) if do_rc else None
        else:
            pos = np.arange(n_pos, dtype=np.intp)
            rc  = (np.int64(4 ** k - 1) - rci) if do_rc else None

        return pos, fwd, rc

    enc_a = enc_table[np.frombuffer(seq_a.encode('ascii'), dtype=np.uint8)]
    enc_b = enc_table[np.frombuffer(seq_b.encode('ascii'), dtype=np.uint8)]
    pos_a, fwd_a, rc_a = _kmer_hashes(enc_a, do_rc=True)
    pos_b, fwd_b, _    = _kmer_hashes(enc_b, do_rc=False)  # only fwd for seq_b

    bin_a = (pos_a // binsize).astype(np.int64)
    bin_b = (pos_b // binsize).astype(np.int64)

    # Pack (hash << SHIFT | bin) into a single int64 so np.sort (no index array)
    # replaces argsort + fancy-index.  Sorting the b-array also makes searchsorted
    # queries sequential → cache-friendly (avoids random jumps through a_packed).
    # SHIFT bits hold the bin index; hash occupies the upper bits.
    # Constraint: 2*k + SHIFT < 63 (int64 must not overflow).
    # For k<=20 SHIFT=21 supports up to 2M bins (~2Gbp at binsize=1000).
    SHIFT = max(int(max(bins_a, bins_b) - 1).bit_length(), 1)
    if 2 * k + SHIFT >= 63:
        # Fall back to argsort for unusually large k or tiny binsize.
        so_a    = np.argsort(np.concatenate([fwd_a, rc_a]), kind='stable')
        a_h_all = np.concatenate([fwd_a, rc_a])[so_a]
        a_b_all = np.concatenate([bin_a, bin_a])[so_a]
        so_b    = np.argsort(fwd_b, kind='stable')
        fwd_b_s = fwd_b[so_b]
        bin_b_s = bin_b[so_b]
        lo = np.searchsorted(a_h_all, fwd_b_s, side='left')
        hi = np.searchsorted(a_h_all, fwd_b_s, side='right')
        MASK = None
        def _get_ba(pk_arr, l, h): return a_b_all[l:h]
        def _get_ba_rep(lo_m, cnt_m):
            grp = np.empty(len(lo_m), dtype=np.int64); grp[0] = 0
            if len(lo_m) > 1: grp[1:] = np.cumsum(cnt_m[:-1])
            diffs = np.ones(int(cnt_m.sum()), dtype=np.int64); diffs[0] = lo_m[0]
            if len(lo_m) > 1: diffs[grp[1:]] = lo_m[1:] - lo_m[:-1] - cnt_m[:-1] + 1
            return a_b_all[np.cumsum(diffs)]
    else:
        MASK  = (1 << SHIFT) - 1
        a_h_c = np.concatenate([fwd_a, rc_a])
        a_b_c = np.concatenate([bin_a, bin_a])
        a_pk  = (a_h_c << SHIFT) | a_b_c;  a_pk.sort()

        # Repeat-filter: drop k-mers whose hash appears > max_kmer_freq times in
        # seq_a (satellite DNA, transposons, …).  These would explode total matches
        # and produce only uninformative noise in the plot.
        if max_kmer_freq > 0 and len(a_pk) > 0:
            a_hash_part = a_pk >> SHIFT
            run_ends    = np.concatenate([np.where(np.diff(a_hash_part))[0] + 1, [len(a_pk)]])
            run_starts  = np.concatenate([[0], run_ends[:-1]])
            run_lens    = run_ends - run_starts
            valid_runs  = run_lens <= max_kmer_freq
            keep        = np.zeros(len(a_pk) + 1, dtype=np.int8)
            keep[run_starts[valid_runs]] += 1
            keep[run_ends[valid_runs]]   -= 1
            a_pk = a_pk[np.cumsum(keep)[:-1].astype(bool)]
            n_filtered = valid_runs.size - valid_runs.sum()
            if n_filtered:
                print(f"  [dotplot_kmer] filtered {n_filtered} high-freq k-mer hashes "
                      f"(>{max_kmer_freq} occurrences)")

        # Collapse duplicate (hash, bin_a) pairs: multiple copies of the same k-mer
        # in the same bin (tandem repeats, satellite DNA) reduce to one entry.
        # a_pk is already sorted so np.unique is a O(n) linear scan.
        # Without this, a k-mer appearing N times in bin_a and M times in bin_b
        # would contribute N*M identical flat indices — pure wasted work.
        a_pk = np.unique(a_pk)

        b_pk  = (fwd_b << SHIFT) | bin_b;   b_pk.sort()
        b_pk  = np.unique(b_pk)
        fwd_b_s = b_pk >> SHIFT
        bin_b_s = b_pk & MASK
        lo = np.searchsorted(a_pk, fwd_b_s << SHIFT,       side='left')
        hi = np.searchsorted(a_pk, (fwd_b_s + 1) << SHIFT, side='left')
        def _get_ba_rep(lo_m, cnt_m):
            grp = np.empty(len(lo_m), dtype=np.int64); grp[0] = 0
            if len(lo_m) > 1: grp[1:] = np.cumsum(cnt_m[:-1])
            total = int(cnt_m.sum())
            diffs = np.ones(total, dtype=np.int64); diffs[0] = lo_m[0]
            if len(lo_m) > 1: diffs[grp[1:]] = lo_m[1:] - lo_m[:-1] - cnt_m[:-1] + 1
            return a_pk[np.cumsum(diffs)] & MASK

    cnt = hi - lo
    has_match = cnt > 0
    mat = np.zeros((bins_a, bins_b), dtype=np.uint32)
    import time as _t
    _tt = _t.perf_counter()
    if has_match.any():
        lo_m  = lo[has_match]
        cnt_m = cnt[has_match]
        bb_m  = bin_b_s[has_match]
        total_matches = int(cnt_m.sum())
        print(f"  [dotplot_kmer] total k-mer matches to accumulate: {total_matches:,}")
        ba_rep = _get_ba_rep(lo_m, cnt_m)
        print(f"  [dotplot_kmer]   _get_ba_rep: {_t.perf_counter()-_tt:.2f}s"); _tt=_t.perf_counter()
        bb_rep = np.repeat(bb_m, cnt_m)
        print(f"  [dotplot_kmer]   np.repeat: {_t.perf_counter()-_tt:.2f}s"); _tt=_t.perf_counter()
        flat   = ba_rep * bins_b + bb_rep
        print(f"  [dotplot_kmer]   flat: {_t.perf_counter()-_tt:.2f}s"); _tt=_t.perf_counter()
        mat    = np.bincount(flat, minlength=bins_a * bins_b).reshape(bins_a, bins_b).astype(np.uint32)
        print(f"  [dotplot_kmer]   bincount: {_t.perf_counter()-_tt:.2f}s"); _tt=_t.perf_counter()

    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))

    plot_mat = np.log1p(mat) if log_scale else mat
    nz = plot_mat[plot_mat > 0]
    vmin = 0.0
    vmax = float(np.percentile(nz, 99.5)) if nz.size > 0 else 1.0
    print(f"  [dotplot_kmer]   vmin/vmax: {_t.perf_counter()-_tt:.2f}s"); _tt=_t.perf_counter()
    ax.imshow(plot_mat, origin='lower', aspect='equal', cmap='Greys', vmin=vmin, vmax=vmax)
    print(f"  [dotplot_kmer]   imshow: {_t.perf_counter()-_tt:.2f}s")
    ax.set_xlabel("Seq B bins")
    ax.set_ylabel("Seq A bins")
    ax.set_title(f"K-mer dotplot (k={k}, binsize={binsize})")
    return ax


def _bins_to_genome_intervals(bin_values: list, segs_info: list, binsize: int) -> list:
    """
    Map hot bins (bin_values[i] > 0) back to merged genome (contig, start, end) intervals.
    segs_info: [(contig, seg_start, seg_end, plot_offset), ...] from _build_axis_map.
    """
    intervals = []
    for i, v in enumerate(bin_values):
        if v <= 0:
            continue
        plot_start = i * binsize
        for contig, seg_start, seg_end, offset in segs_info:
            seg_plot_end = offset + math.ceil((seg_end - seg_start) / binsize) * binsize
            if offset <= plot_start < seg_plot_end:
                gen_start = seg_start + (plot_start - offset)
                gen_end   = min(seg_end, seg_start + ((i + 1) * binsize - offset))
                intervals.append((contig, gen_start, gen_end))
                break

    # Merge consecutive intervals on the same contig
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for contig, s, e in intervals[1:]:
        if contig == merged[-1][0] and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([contig, s, e])
    return [(c, s, e) for c, s, e in merged]


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


def reads_alignment_dotplot(sc: SegmentCollection, read_names: set[str], pid_threshold: float,
                            plen_threshold: float, db_path: str, ax=None, line_alpha: float = 1,
                            line_width_px: float = 6.0, pid_tie_eps: float = PID_TIE_EPS_DEFAULT,
                            boundary_color: str = "black", boundary_alpha: float = 1.0,
                            boundary_lw: float = 1.0, asm1: str = "asm1", asm2: str = "asm2",
                            aln_by_asm: dict | None = None, binsize: int = 1000,
                            read_classes: dict | None = None):
    """
    Plot read support as a dot-plot (lines) between asm2 (X) and asm1 (Y).
    Alignments are fetched from the read-indexed alignment database.

    If aln_by_asm is provided it is used directly and no DB query is issued.
    """
    def _map_pos(lookup: dict, contig: str, start: int, end: int):
        for s, e, off in lookup.get(contig, ()):
            if start >= s and end <= e:
                return off + (start - s), off + (end - s)
        return None

    standalone = ax is None
    if standalone:
        _, ax = plt.subplots(figsize=(16, 5))

    sc = sorted(sc, key=lambda s: (s.contig, s.start, s.end))
    asm1_segs, asm1_len = _build_axis_map(sc, asm1, binsize=binsize)
    asm2_segs, asm2_len = _build_axis_map(sc, asm2, binsize=binsize)

    # Build contig-keyed lookups for O(1) contig access in _map_pos.
    asm1_lookup: dict = {}
    for contig, s, e, off in asm1_segs:
        asm1_lookup.setdefault(contig, []).append((s, e, off))
    asm2_lookup: dict = {}
    for contig, s, e, off in asm2_segs:
        asm2_lookup.setdefault(contig, []).append((s, e, off))

    if not asm1_segs or not asm2_segs:
        raise ValueError("SegmentCollection must contain both asm1 and asm2 segments")

    reads = list(read_names)
    if not reads:
        return ax

    if aln_by_asm is None:
        con = sqlite3.connect(db_path)
        try:
            aln_by_asm = _fetch_aln_from_db(con, reads, [asm1, asm2], pid_threshold, plen_threshold)
        finally:
            con.close()

    # Group alignments by read.
    by_read = {}
    for asm, alns in aln_by_asm.items():
        for read_name, contig, start, end, strand, pid, plen in alns:
            by_read.setdefault(read_name, {}).setdefault(asm, []).append(
                (contig, start, end, strand, pid, plen)
            )

    # Build one plotted line per read using best-PID alignment on each assembly.
    lines = []
    for read_name in reads:
        per_asm = by_read.get(read_name, {})
        asm1_alns = per_asm.get(asm1, [])
        asm2_alns = per_asm.get(asm2, [])
        if not asm1_alns or not asm2_alns:
            continue
        r_contig, r_start, r_end, r_strand, r_pid, r_plen = max(
            asm1_alns, key=lambda x: (x[4], x[5], x[2] - x[1])
        )
        c_contig, c_start, c_end, c_strand, c_pid, c_plen = max(
            asm2_alns, key=lambda x: (x[4], x[5], x[2] - x[1])
        )

        r_mapped = _map_pos(asm1_lookup, r_contig, r_start, r_end)
        c_mapped = _map_pos(asm2_lookup, c_contig, c_start, c_end)
        if r_mapped is None or c_mapped is None:
            continue

        r0, r1 = r_mapped
        c0, c1 = c_mapped
        if r_strand == 0:
            r0, r1 = r1, r0
        if c_strand == 0:
            c0, c1 = c1, c0
        lines.append((c0, c1, r0, r1, c_pid, r_pid, read_name))

    # Plot fixed-width lines with PID-based colors, batched by color into LineCollections.
    if lines:
        dpi = ax.figure.dpi if ax.figure is not None else 72.0
        lw_points = line_width_px * 72.0 / dpi
        asm1_segs_lc, asm2_segs_lc = [], []
        for c0, c1, r0, r1, c_pid, r_pid, rn in lines:
            if read_classes is not None:
                if rn in read_classes.get("asm1_best", ()):
                    asm1_segs_lc.append([(c0, r0), (c1, r1)])
                elif rn in read_classes.get("asm2_best", ()):
                    asm2_segs_lc.append([(c0, r0), (c1, r1)])
                # asm1_only/asm2_only only map to one side — no line possible
                # common/other — omit from overlay (no informative direction)
            elif c_pid - r_pid > pid_tie_eps:
                asm2_segs_lc.append([(c0, r0), (c1, r1)])
            elif r_pid - c_pid > pid_tie_eps:
                asm1_segs_lc.append([(c0, r0), (c1, r1)])
        for segs_lc, color in (
            (asm1_segs_lc, SUPPORT_COLORS["asm1_best"]),
            (asm2_segs_lc, SUPPORT_COLORS["asm2_best"]),
        ):
            if segs_lc:
                ax.add_collection(LineCollection(
                    segs_lc, colors=color, linewidths=lw_points, alpha=line_alpha
                ))

    # Segment boundaries: with bin-aligned offsets each segment starts at an
    # exact multiple of binsize, so the boundary between seg j and seg j+1 is
    # simply seg[j+1].offset — guaranteed on a bin edge, no snapping needed.
    # Same-contig boundaries are drawn much paler than cross-contig boundaries.
    for i in range(1, len(asm2_segs)):
        same_contig = asm2_segs[i][0] == asm2_segs[i - 1][0]
        ax.axvline(asm2_segs[i][3],
                   color=boundary_color,
                   alpha=boundary_alpha * (0.5 if same_contig else 1.0),
                   linewidth=boundary_lw * (0.5 if same_contig else 1.0))
    for i in range(1, len(asm1_segs)):
        same_contig = asm1_segs[i][0] == asm1_segs[i - 1][0]
        ax.axhline(asm1_segs[i][3],
                   color=boundary_color,
                   alpha=boundary_alpha * (0.5 if same_contig else 1.0),
                   linewidth=boundary_lw * (0.5 if same_contig else 1.0))

    # Segment captions (major ticks, no tick lines) — coordinates only
    x_centers = [off + (end - start) / 2 for _, start, end, off in asm2_segs]
    x_labels = [f"{start:,}-{end:,}" for _, start, end, _ in asm2_segs]
    y_centers = [off + (end - start) / 2 for _, start, end, off in asm1_segs]
    y_labels = [f"{start:,}-{end:,}" for _, start, end, _ in asm1_segs]
    ax.set_xticks(x_centers)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_yticks(y_centers)
    ax.set_yticklabels(y_labels, rotation=0, ha="right", va="center")

    # Coordinate ticks every 50,000 bp (minor ticks)
    x_minor = np.arange(0, asm2_len + 1, 50_000)
    y_minor = np.arange(0, asm1_len + 1, 50_000)
    ax.xaxis.set_minor_locator(mticker.FixedLocator(x_minor))
    ax.yaxis.set_minor_locator(mticker.FixedLocator(y_minor))
    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

    ax.tick_params(axis="x", which="major", length=0, pad=20)
    ax.tick_params(axis="y", which="major", length=0, pad=6)
    ax.tick_params(axis="x", which="minor", length=3, width=2, labelsize=8, labelbottom=True)
    ax.tick_params(axis="y", which="minor", length=3, width=2, labelleft=False)

    ax.text(1.0, -0.08, "Mbp", transform=ax.transAxes, ha="right", va="top")
    ax.set_xlabel(f"{asm2} (concatenated segments)")
    ax.set_ylabel(f"{asm1} (concatenated segments)")

    if standalone:
        ax.set_xlim(0, max(asm2_len - 1, 0))
        ax.set_ylim(0, max(asm1_len - 1, 0))
        ax.set_aspect('equal', adjustable='box')
    return ax


def dotplot_kmer_mosaic(seqs_a: List[str], seqs_b: List[str],
                        k: int = 15, binsize: int = 1000, ax=None,
                        max_cells: int = 10_000_000,
                        mask_ambiguous: bool = True, log_scale: bool = False,
                        max_kmer_freq: int = 100):
    """
    K-mer dotplot rendered as a mosaic of independent per-segment-pair sub-dotplots.

    For each pair (seqs_a[i], seqs_b[j]) a separate k-mer match matrix is computed
    and placed in the corresponding cell of the combined grid.  Because every segment
    starts on a bin boundary (offset = cumulative ceil(seg_len/binsize)*binsize), no
    bin ever straddles two segments — eliminating the "mixed bin" artefact.

    Coordinate space: x ∈ [0, Σ_j ceil(len_j/binsize)*binsize],
                      y ∈ [0, Σ_i ceil(len_i/binsize)*binsize].
    Each pixel is exactly binsize × binsize bp.
    """
    if k <= 0 or binsize <= 0:
        raise ValueError("k and binsize must be positive")
    if not seqs_a or not seqs_b:
        raise ValueError("seqs_a and seqs_b must be non-empty")

    enc_table = np.full(256, 255, dtype=np.uint8)
    for _b, _v in zip(b'ACGT', [0, 1, 2, 3]):
        enc_table[_b] = _v

    bins_a = [max(math.ceil(len(s) / binsize), 1) for s in seqs_a]
    bins_b = [max(math.ceil(len(s) / binsize), 1) for s in seqs_b]
    total_bins_a = sum(bins_a)
    total_bins_b = sum(bins_b)

    if total_bins_a * total_bins_b > max_cells:
        raise ValueError(
            f"Mosaic matrix {total_bins_a}×{total_bins_b} = {total_bins_a*total_bins_b} cells; "
            "increase binsize or max_cells"
        )

    max_bins_any = max(max(bins_a), max(bins_b))
    SHIFT = max(int(max_bins_any - 1).bit_length(), 1)
    if 2 * k + SHIFT >= 63:
        raise ValueError(
            f"k={k} too large for binsize={binsize} (SHIFT={SHIFT}); reduce k or increase binsize"
        )
    MASK = (1 << SHIFT) - 1

    def _seg_pk(seq: str, do_rc: bool) -> np.ndarray:
        """Packed sorted (hash << SHIFT | bin) array for one segment sequence."""
        seq = seq.upper()
        if len(seq) < k:
            return np.array([], dtype=np.int64)
        enc = enc_table[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]
        n_pos = len(enc) - k + 1
        ph = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
        pl = 4 ** np.arange(k,           dtype=np.int64) if do_rc else None
        ei = enc.astype(np.int64)
        fwd = np.zeros(n_pos, dtype=np.int64)
        rci = np.zeros(n_pos, dtype=np.int64) if do_rc else None
        for j in range(k):
            sl = ei[j: j + n_pos]
            fwd += sl * ph[j]
            if do_rc:
                rci += sl * pl[j]
        if mask_ambiguous:
            inv_pos = np.where(enc == 255)[0]
            valid = np.ones(n_pos, dtype=bool)
            for p in inv_pos:
                valid[max(0, p - k + 1): min(n_pos, p + 1)] = False
            pos = np.where(valid)[0]
            fwd = fwd[valid]
            rc  = (np.int64(4 ** k - 1) - rci[valid]) if do_rc else None
        else:
            pos = np.arange(n_pos, dtype=np.intp)
            rc  = (np.int64(4 ** k - 1) - rci) if do_rc else None
        seg_bins = (pos // binsize).astype(np.int64)
        if do_rc:
            h_all = np.concatenate([fwd, rc])
            b_all = np.concatenate([seg_bins, seg_bins])
        else:
            h_all = fwd
            b_all = seg_bins
        pk = (h_all << SHIFT) | b_all
        pk.sort()
        # Repeat filter on a-side
        if do_rc and max_kmer_freq > 0 and len(pk) > 0:
            hp = pk >> SHIFT
            re = np.concatenate([np.where(np.diff(hp))[0] + 1, [len(pk)]])
            rs = np.concatenate([[0], re[:-1]])
            rl = re - rs
            vr = rl <= max_kmer_freq
            keep = np.zeros(len(pk) + 1, dtype=np.int8)
            keep[rs[vr]] += 1
            keep[re[vr]] -= 1
            pk = pk[np.cumsum(keep)[:-1].astype(bool)]
        return np.unique(pk)   # collapse duplicate (hash, bin) pairs

    # Hash each segment once
    pks_a = [_seg_pk(s, do_rc=True)  for s in seqs_a]   # fwd + RC
    pks_b = [_seg_pk(s, do_rc=False) for s in seqs_b]   # fwd only

    cum_a = [0] + list(np.cumsum(bins_a))
    cum_b = [0] + list(np.cumsum(bins_b))

    mat = np.zeros((total_bins_a, total_bins_b), dtype=np.uint32)

    for i, a_pk in enumerate(pks_a):
        if len(a_pk) == 0:
            continue
        for j, b_pk in enumerate(pks_b):
            if len(b_pk) == 0:
                continue
            b_fwd_s = b_pk >> SHIFT
            bin_b_s = b_pk & MASK
            lo = np.searchsorted(a_pk, b_fwd_s << SHIFT,       side='left')
            hi = np.searchsorted(a_pk, (b_fwd_s + 1) << SHIFT, side='left')
            cnt = hi - lo
            has_match = cnt > 0
            if not has_match.any():
                continue
            lo_m  = lo[has_match]
            cnt_m = cnt[has_match]
            bb_m  = bin_b_s[has_match]
            total = int(cnt_m.sum())
            # Expand a_pk bin indices (ragged-range cumsum)
            grp = np.empty(len(lo_m), dtype=np.int64); grp[0] = 0
            if len(lo_m) > 1: grp[1:] = np.cumsum(cnt_m[:-1])
            diffs = np.ones(total, dtype=np.int64); diffs[0] = lo_m[0]
            if len(lo_m) > 1:
                diffs[grp[1:]] = lo_m[1:] - lo_m[:-1] - cnt_m[:-1] + 1
            ba_rep = a_pk[np.cumsum(diffs)] & MASK
            bb_rep = np.repeat(bb_m, cnt_m)
            flat   = ba_rep * bins_b[j] + bb_rep
            sub    = np.bincount(flat, minlength=bins_a[i] * bins_b[j])
            mat[cum_a[i]:cum_a[i+1], cum_b[j]:cum_b[j+1]] = \
                sub.reshape(bins_a[i], bins_b[j]).astype(np.uint32)

    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))

    plot_mat = np.log1p(mat) if log_scale else mat.astype(float)
    nz = plot_mat[plot_mat > 0]
    vmin = 0.0
    vmax = float(np.percentile(nz, 99.5)) if nz.size > 0 else 1.0
    ax.imshow(plot_mat, origin='lower', aspect='equal', cmap='Greys',
              vmin=vmin, vmax=vmax,
              extent=(0, total_bins_b * binsize, 0, total_bins_a * binsize))
    ax.set_xlabel("Seq B bins")
    ax.set_ylabel("Seq A bins")
    ax.set_title(f"K-mer dotplot mosaic (k={k}, binsize={binsize})")
    return ax


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


def likelihood_report(db_path: str, read_classes: dict,
                      asm1_segs: list, asm2_segs: list,
                      asm1: str, asm2: str,
                      pid_threshold: float = 0.0,
                      plen_threshold: float = 0.0) -> None:
    """
    Compute and print read-support likelihood statistics.

    Pr(r | A) = (1 / L_A) * Σ_p exp(-D(r, p_i))
    where D is the edit distance (NM tag) and L_A is the total length of
    recruited segments for assembly A.

    Likelihood is computed over reads in: asm1_best, asm2_best, common, other.
    asm1_only / asm2_only are counted and reported but excluded from the sum.
    """
    _CLASSES = ["asm1_only", "asm2_only", "asm1_best", "asm2_best", "common", "other"]
    counts = {k: len(read_classes.get(k, set())) for k in _CLASSES}

    L_A1 = sum(s.end - s.start for s in asm1_segs)
    L_A2 = sum(s.end - s.start for s in asm2_segs)

    included = (read_classes.get("asm1_best", set()) |
                read_classes.get("asm2_best", set()) |
                read_classes.get("common",    set()) |
                read_classes.get("other",     set()))

    plen_thr = plen_threshold if plen_threshold is not None else 0.0

    ed_by_asm: dict = {asm1: {}, asm2: {}}
    if included:
        reads_list = list(included)
        con = sqlite3.connect(db_path)
        try:
            for i in range(0, len(reads_list), 800):
                chunk = reads_list[i:i + 800]
                reads_ph = ",".join(["?"] * len(chunk))
                sql = f"""
                    SELECT assembly, read_name, ed FROM aln
                    WHERE assembly IN (?, ?)
                      AND contig IS NOT NULL
                      AND pid >= ? AND plen >= ?
                      AND ed IS NOT NULL
                      AND read_name IN ({reads_ph})
                """
                for row_asm, rname, ed in con.execute(
                        sql, [asm1, asm2, pid_threshold, plen_thr] + chunk):
                    if row_asm in ed_by_asm:
                        ed_by_asm[row_asm].setdefault(rname, []).append(int(ed))
        finally:
            con.close()

    def _logsumexp_neg(eds):
        m = min(eds)
        return -m + math.log(sum(math.exp(-(e - m)) for e in eds))

    def _total_log_likelihood(asm_ed_map, L_A):
        if L_A <= 0:
            return float('-inf'), len(included)
        ll, n_missing = 0.0, 0
        log_L_A = math.log(L_A)
        for rname in included:
            eds = asm_ed_map.get(rname)
            if not eds:
                n_missing += 1
                continue
            ll += _logsumexp_neg(eds) - log_L_A
        return ll, n_missing

    ll1, miss1 = _total_log_likelihood(ed_by_asm[asm1], L_A1)
    ll2, miss2 = _total_log_likelihood(ed_by_asm[asm2], L_A2)
    n_inc = len(included)
    delta = ll2 - ll1 if math.isfinite(ll1) and math.isfinite(ll2) else float('nan')

    w = max(len(asm1), len(asm2), 12)
    print()
    print("=== Likelihood Report ===")
    print("Read class counts:")
    for cls in _CLASSES:
        print(f"  {cls:<12} : {counts[cls]:>7,}")
    print(f"  {'total':<12} : {sum(counts.values()):>7,}")
    print()
    print("Recruited segment lengths:")
    print(f"  {asm1:<{w}} : {L_A1:>14,} bp")
    print(f"  {asm2:<{w}} : {L_A2:>14,} bp")
    print()
    def _fmt_coords(coord_map):
        parts = [f"{c}:{s}-{e}({n})" for c, ivs in coord_map.items() for s, e, n in ivs]
        return "[" + ", ".join(parts) + "]" if parts else "[]"

    asm1_coords = _stripe_coords(db_path, read_classes.get("asm1_only", set()), asm1, pid_threshold, plen_threshold)
    asm2_coords = _stripe_coords(db_path, read_classes.get("asm2_only", set()), asm2, pid_threshold, plen_threshold)
    print(f"{asm1} support only: {_fmt_coords(asm1_coords)}")
    print(f"{asm2} support only: {_fmt_coords(asm2_coords)}")
    print()
    print(f"Log-likelihood  (n={n_inc:,} reads: asm1_best + asm2_best + common + other)")
    print(f"  log P(reads | {asm1:<{w}}) = {ll1:>14.3f}   [reads without aln: {miss1}]")
    print(f"  log P(reads | {asm2:<{w}}) = {ll2:>14.3f}   [reads without aln: {miss2}]")
    if math.isfinite(delta):
        better = asm2 if delta > 0 else asm1
        sign = "+" if delta > 0 else ""
        print(f"  ΔLLR  ({asm2} − {asm1}){' ' * (w - len(asm2))} = {sign}{delta:>13.3f}   [{better} better supported]")
    else:
        print(f"  ΔLLR  ({asm2} − {asm1}) = N/A")
    print()


def _coord_ticks(segs_info, total_len):
    """Major and minor tick positions for genome coordinates across concatenated segments.
    Major ticks: one per segment at its start (segment boundary).
    Minor ticks: absolute genome multiples of step inside each segment; step auto-selected
    to give ≤20 ticks across total_len.
    Returns (major_positions, minor_positions) where minor excludes major positions."""
    _STEPS = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000,
              100_000, 200_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000]
    minor_step = _STEPS[-1]
    for s in _STEPS:
        if total_len / s <= 40:
            minor_step = s
            break

    major_pos = [float(off) for _, _, _, off in segs_info]
    major_set = set(major_pos)

    minor_pos = []
    for _, gen_s, gen_e, off in segs_info:
        first = math.ceil(gen_s / minor_step) * minor_step
        for g in range(first, gen_e, minor_step):
            plot_coord = off + (g - gen_s)
            if 0 < plot_coord <= total_len and float(plot_coord) not in major_set:
                minor_pos.append(float(plot_coord))

    return major_pos, minor_pos


def report_plot(seq_asm1: str, seq_asm2: str, segments: SegmentCollection, read_names: set[str],
                db_path: str, pid_threshold: float, plen_threshold: float, k: int = 15, binsize: int | None = None,
                output_path: str = "report_plot.png", dpi: int = 150,
                asm2_bam: str | None = None, asm1_bam: str | None = None,
                asm1_fasta: str | None = None, asm2_fasta: str | None = None,
                boundary_params: tuple[str, float, float] = ("black", 1.0, 1.0),
                asm1: str = "asm1", asm2: str = "asm2",
                pid_tie_eps: float = PID_TIE_EPS_DEFAULT,
                best_pid_delta: float = BEST_PID_DELTA_DEFAULT,
                max_kmer_freq: int = 100,
                max_coverage: float | None = None, show: bool = True,
                support_stripes: bool = True,
                seq_segments: SegmentCollection | None = None,
                title: str | None = None,
                classification: str = "support"):
    """
    Create a combined report plot with k-mer dotplot as background and
    read alignment dotplot overlaid on top. Saves to output_path.

    Parameters (marginal plots):
    - max_coverage: fix the coverage axis limit for both marginal alignment plots
      (ax_top y-axis and ax_right x-axis) to this value.  When None (default),
      the limit is set to max(maxx, maxy) derived from the data, which makes both
      marginals share the same scale within a single call but may differ across
      different calls.  Pass the same value to all report_plot calls to make
      plots directly comparable.

    - show: when True, display the figure with plt.show(); set False for batch runs.
    """
    import time as _time
    def _t(label, t0):
        print(f"  [report_plot] {label}: {_time.perf_counter()-t0:.2f}s")
        return _time.perf_counter()

    _t0 = _time.perf_counter()

    sc_sorted = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    # seq_sc defines the coordinate system for the sequences (seq_asm1/seq_asm2),
    # the dotplot extent, alignment line mapping, and contig bars.
    # When seq_segments is provided (sequences were built from different segments than
    # the custom segments passed as `segments`), use it here to keep alignment lines
    # consistent with the dotplot content.  `segments` is kept only for marginal
    # BAM coverage queries.
    seq_sc = (sorted(seq_segments, key=lambda s: (s.contig, s.start, s.end))
              if seq_segments is not None else sc_sorted)

    # Auto-compute binsize from segment lengths when not provided:
    # target ~1000 bins on the longer axis, but never smaller than k.
    if binsize is None:
        total_a = sum(seg.end - seg.start for seg in seq_sc if seg.assembly == asm1)
        total_b = sum(seg.end - seg.start for seg in seq_sc if seg.assembly == asm2)
        binsize = max(k, math.ceil(max(total_a, total_b, 1) / 1000))

    # Bin-aligned axis maps: each segment padded to a whole number of bins so
    # offsets (and therefore segment boundaries) are exact multiples of binsize.
    asm1_segs_info, plot_h = _build_axis_map(seq_sc, asm1, binsize=binsize)
    asm2_segs_info, plot_w = _build_axis_map(seq_sc, asm2, binsize=binsize)

    # Segment lists used by marginal plots and BAM read collection (always from segments).
    asm2_segs_list = [s for s in sc_sorted if s.assembly == asm2]
    asm1_segs_list = [s for s in sc_sorted if s.assembly == asm1]

    # Fetch sequences internally from seq_sc order when fasta paths are provided.
    # This guarantees the sequence order always matches the axis map, regardless of
    # how the caller built any pre-passed seq_asm1/seq_asm2 strings.
    if asm1_fasta is not None and asm2_fasta is not None:
        _fa1 = pysam.FastaFile(asm1_fasta)
        _fa2 = pysam.FastaFile(asm2_fasta)
        try:
            _p1, _p2 = [], []
            for seg in seq_sc:
                if seg.assembly == asm1:
                    _p1.append(_fa1.fetch(seg.contig, seg.start, seg.end))
                elif seg.assembly == asm2:
                    _p2.append(_fa2.fetch(seg.contig, seg.start, seg.end))
            seq_asm1 = "".join(_p1)
            seq_asm2 = "".join(_p2)
        finally:
            _fa1.close()
            _fa2.close()

    _t0 = _t(f"setup  seqs={len(seq_asm1)//1000}kbp/{len(seq_asm2)//1000}kbp  reads={len(read_names)}  segs={len(sc_sorted)}", _t0)

    # Expand read set: add all BAM reads that touch the segments so that
    # classification covers every read visible in the marginal coverage plots.
    all_reads: set = set(read_names)
    if asm2_bam is not None:
        all_reads |= _collect_bam_read_names(asm2_bam, asm2_segs_list, pid_threshold, plen_threshold)
    if asm1_bam is not None:
        all_reads |= _collect_bam_read_names(asm1_bam, asm1_segs_list, pid_threshold, plen_threshold)
    _t0 = _t(f"collect_bam_read_names  ({len(all_reads) - len(read_names):+d} reads)", _t0)

    # Fetch alignments once for all reads; reused by dotplot overlay and classifier.
    con = sqlite3.connect(db_path)
    try:
        aln_by_asm = _fetch_aln_from_db(con, list(all_reads), [asm1, asm2], pid_threshold, plen_threshold)
    finally:
        con.close()
    _t0 = _t(f"fetch_aln_from_db  ({sum(len(v) for v in aln_by_asm.values())} rows)", _t0)

    # Classify reads early so the overlay and marginals share the same read_classes.
    # segs_by_asm=None → use all alignments globally so that reads shown in a
    # subset view keep their globally-correct category (avoids misclassification
    # when a read's best asm1 alignment falls outside the displayed subset).
    read_classes = None
    if asm2_bam is not None and asm1_bam is not None:
        fa1 = pysam.FastaFile(asm1_fasta) if asm1_fasta else None
        fa2 = pysam.FastaFile(asm2_fasta) if asm2_fasta else None
        try:
            read_classes = _classify_reads_for_support_colors(
                db_path, all_reads, pid_threshold, plen_threshold,
                segs_by_asm=None, asm1=asm1, asm2=asm2, pid_tie_eps=pid_tie_eps,
                aln_by_asm=aln_by_asm, asm1_fasta=fa1, asm2_fasta=fa2,
                best_pid_delta=best_pid_delta,
            )
        finally:
            if fa1: fa1.close()
            if fa2: fa2.close()
        _t0 = _t("classify_reads_for_support_colors", _t0)

    # Split concatenated sequences into per-segment lists for the mosaic dotplot.
    # The split uses raw segment lengths (not bin-padded), matching the order of
    # asm*_segs_info which was built from sc_sorted.
    seqs_a, raw_off = [], 0
    for _, start, end, _ in asm1_segs_info:
        seqs_a.append(seq_asm1[raw_off: raw_off + (end - start)])
        raw_off += end - start
    seqs_b, raw_off = [], 0
    for _, start, end, _ in asm2_segs_info:
        seqs_b.append(seq_asm2[raw_off: raw_off + (end - start)])
        raw_off += end - start

    # Figure height: main axes is bins_a/bins_b tall for square dotplot pixels,
    # +2 inches for ax_top marginal plot.
    total_bins_a = plot_h // binsize
    total_bins_b = plot_w // binsize
    height = max(4.0, 16.0 * total_bins_a / max(total_bins_b, 1) + 2)
    fig, ax = plt.subplots(figsize=(16.0, height), dpi=dpi)

    # Background k-mer mosaic dotplot (no mixed bins at segment boundaries)
    dotplot_kmer_mosaic(seqs_a, seqs_b, k=k, binsize=binsize, ax=ax,
                        max_kmer_freq=max_kmer_freq)
    ax.set_aspect('auto')
    _t0 = _t("dotplot_kmer_mosaic", _t0)

    # Full-size inset for reads overlay
    ax_overlay = ax.inset_axes((0, 0, 1, 1))
    ax_overlay.set_facecolor("none")
    ax_overlay.patch.set_alpha(0.0)
    ax_overlay.set_frame_on(False)
    b_color, b_alpha, b_lw = boundary_params
    reads_alignment_dotplot(
        seq_sc, all_reads, pid_threshold, plen_threshold, db_path,
        ax=ax_overlay, boundary_color=b_color, boundary_alpha=b_alpha, boundary_lw=b_lw,
        asm1=asm1, asm2=asm2, pid_tie_eps=pid_tie_eps, aln_by_asm=aln_by_asm,
        binsize=binsize, read_classes=read_classes,
    )
    _t0 = _t("reads_alignment_dotplot", _t0)

    # Keep overlays aligned with base axes
    ax.set_xlim(0, plot_w)
    ax.set_ylim(0, plot_h)
    ax_overlay.set_xlim(ax.get_xlim())
    ax_overlay.set_ylim(ax.get_ylim())

    # Segment color bars — one colored rectangle per segment, colored by contig.
    # Strips sit between the main dotplot and the marginal alignment plots.
    all_contigs = sorted({c for c, _, _, _ in asm1_segs_info + asm2_segs_info})
    _cmap = plt.cm.get_cmap('tab10' if len(all_contigs) <= 10 else 'tab20')
    contig_color = {c: _cmap(i % _cmap.N) for i, c in enumerate(all_contigs)}

    ax_top = ax_right = None   # set inside the marginal-plots block if bams are provided

    divider = make_axes_locatable(ax)
    ax_xbar = divider.append_axes("bottom", size=0.15, pad=0.0, sharex=ax)
    ax_ybar = divider.append_axes("left",   size=0.15, pad=0.0, sharey=ax)

    # Hide everything on both bar axes except the outer edge that will carry tick marks.
    # axis("off") would suppress the tick marks we want, so hide selectively.
    ax_xbar.set_facecolor("none")
    for sp in ("top", "left", "right"): ax_xbar.spines[sp].set_visible(False)
    ax_xbar.yaxis.set_visible(False)
    ax_xbar.xaxis.set_ticks_position('bottom')  # force outermost edge; sharex default is 'top'
    ax_xbar.tick_params(axis="x", bottom=True, top=False, labelbottom=False, length=4)

    ax_ybar.set_facecolor("none")
    for sp in ("top", "right", "bottom"): ax_ybar.spines[sp].set_visible(False)
    ax_ybar.xaxis.set_visible(False)
    ax_ybar.tick_params(axis="y", left=True, right=False, labelleft=False, length=4)

    # Caption fitting: estimate pixels per data-unit to decide if text fits in bar.
    # Using approximate axis sizes; tight-layout may shift these slightly.
    _bar_fontsize = 7
    _char_px = _bar_fontsize / 72 * dpi * 0.62   # ~pixels per character at this dpi
    _fig_w_px  = 16.0 * dpi                       # approximate full figure width
    _fig_h_px  = (height - 2.15) * dpi            # approximate main-axes height (minus bar+top)
    _legend_entries = {}   # contig -> color, for bars too small to label in-place

    def _contig_runs(segs_info, total):
        """Yield (contig, bar_off, bar_end) merging consecutive same-contig segments."""
        if not segs_info:
            return
        cur_contig, _, _, cur_off = segs_info[0]
        for i in range(1, len(segs_info)):
            if segs_info[i][0] != cur_contig:
                yield cur_contig, cur_off, segs_info[i][3]
                cur_contig, _, _, cur_off = segs_info[i]
        yield cur_contig, cur_off, total

    ax_xbar.set_xlim(0, plot_w)
    ax_xbar.set_ylim(0, 1)
    for contig, bar_off, bar_end in _contig_runs(asm2_segs_info, plot_w):
        bar_w = bar_end - bar_off
        ax_xbar.add_patch(plt.Rectangle(
            (bar_off, 0), bar_w, 1,
            color=contig_color[contig], transform=ax_xbar.transData, clip_on=True
        ))
        if bar_w / plot_w * _fig_w_px >= len(contig) * _char_px * 1.4:
            ax_xbar.text(bar_off + bar_w / 2, 0.5, contig,
                         ha='center', va='center', fontsize=_bar_fontsize,
                         color='black', clip_on=True)
        else:
            _legend_entries[contig] = contig_color[contig]

    ax_ybar.set_xlim(0, 1)
    ax_ybar.set_ylim(0, plot_h)
    for contig, bar_off, bar_end in _contig_runs(asm1_segs_info, plot_h):
        bar_h = bar_end - bar_off
        ax_ybar.add_patch(plt.Rectangle(
            (0, bar_off), 1, bar_h,
            color=contig_color[contig], transform=ax_ybar.transData, clip_on=True
        ))
        if bar_h / plot_h * _fig_h_px >= len(contig) * _char_px * 1.4:
            ax_ybar.text(0.5, bar_off + bar_h / 2, contig,
                         ha='center', va='center', fontsize=_bar_fontsize,
                         rotation=90, color='black', clip_on=True)
        else:
            _legend_entries[contig] = contig_color[contig]

    if _legend_entries:
        from matplotlib.patches import Patch
        _handles = [Patch(color=c, label=n) for n, c in _legend_entries.items()]
        fig.legend(handles=_handles, loc='lower center',
                   ncol=min(len(_legend_entries), 6), fontsize=8,
                   bbox_to_anchor=(0.5, 0), frameon=True)

    # Marginal alignment plots
    if asm2_bam is not None and asm1_bam is not None:
        ax_top = divider.append_axes("top", size=2, pad=0, sharex=ax)
        ax_right = divider.append_axes("right", size=2, pad=0, sharey=ax)
        ax_top.tick_params(axis="x", labelbottom=False)
        ax_right.tick_params(axis="y", labelleft=False)

        likelihood_report(db_path, read_classes, asm1_segs_list, asm2_segs_list,
                          asm1, asm2, pid_threshold, plen_threshold)

        if classification == "occurrence":
            _common = read_classes.get("common", set())
            asm2_buckets = _classify_reads_for_occurrence_colors(aln_by_asm, asm2, _common, best_pid_delta=best_pid_delta)
            asm2_colors  = OCCURRENCE_COLORS
            asm1_buckets = _classify_reads_for_occurrence_colors(aln_by_asm, asm1, _common, best_pid_delta=best_pid_delta)
            asm1_colors  = OCCURRENCE_COLORS
        else:  # "support"
            _bucket_keys   = ["asm1_only", "asm2_only", "asm1_best", "asm2_best", "common", "other"]
            _bucket_colors = [SUPPORT_COLORS[k] for k in _bucket_keys]
            asm2_buckets   = [read_classes[k] for k in _bucket_keys]
            asm2_colors    = _bucket_colors
            asm1_buckets   = [read_classes[k] for k in _bucket_keys]
            asm1_colors    = _bucket_colors

        ax_top, maxx, asm2_bucket_bins = alignment_plot(asm2_bam, asm2_segs_list, pid_threshold=pid_threshold,
                                                         plen_threshold=plen_threshold, binsize=binsize,
                                                         ax=ax_top, read_buckets=asm2_buckets,
                                                         read_bucket_colors=asm2_colors,
                                                         orientation="x", return_max=True,
                                                         best_pid_delta=best_pid_delta)
        _t0 = _t(f"alignment_plot asm2  ({len(asm2_segs_list)} segs)", _t0)

        ax_right, maxy, asm1_bucket_bins = alignment_plot(asm1_bam, asm1_segs_list, pid_threshold=pid_threshold,
                                                           plen_threshold=plen_threshold, binsize=binsize,
                                                           ax=ax_right, read_buckets=asm1_buckets,
                                                           read_bucket_colors=asm1_colors,
                                                           orientation="y", return_max=True,
                                                           best_pid_delta=best_pid_delta)
        _t0 = _t(f"alignment_plot asm1  ({len(asm1_segs_list)} segs)", _t0)

        # Transparent stripes on main dotplot: red vertical = asm2_only bins,
        # blue horizontal = asm1_only bins. bucket index 0 = *_only in both bucket lists.
        # Explicit rectangle geometry (BrokenBarHCollection / PatchCollection) avoids
        # fill_between step/where alignment artefacts.
        # Stripes use support-classification bucket indices and are skipped in occurrence mode.
        if support_stripes and classification == "support":
            # bucket order: asm1_only=0, asm2_only=1, asm1_best=2, asm2_best=3, common=4, other=5
            if asm2_bucket_bins is not None:
                patches = [Rectangle((i * binsize, 0), binsize, plot_h)
                           for i, v in enumerate(asm2_bucket_bins[1]) if v > 0]
                if patches:
                    ax.add_collection(PatchCollection(
                        patches, facecolor=SUPPORT_COLORS["asm2_only"], alpha=0.3, linewidth=0,
                    ))
            if asm1_bucket_bins is not None:
                patches = [Rectangle((0, i * binsize), plot_w, binsize)
                           for i, v in enumerate(asm1_bucket_bins[0]) if v > 0]
                if patches:
                    ax.add_collection(PatchCollection(
                        patches, facecolor=SUPPORT_COLORS["asm1_only"], alpha=0.3, linewidth=0,
                    ))
            _t0 = _t("stripes", _t0)

        m = max_coverage if max_coverage is not None else max(maxx, maxy)
        ax_top.set_ylim(0, m)
        ax_right.set_xlim(0, m)
        # Compute shared tick positions so both coverage axes are identical.
        cov_ticks = mticker.MaxNLocator(nbins=5, integer=True).tick_values(0, m)
        cov_ticks = [t for t in cov_ticks if 0 <= t <= m]
        ax_top.set_yticks(cov_ticks)
        ax_right.set_xticks(cov_ticks)
        ax_top.set_xticks([])
        ax_top.set_xlabel("")
        ax_right.set_yticks([])
        ax_right.set_ylabel("")

    # Clear ticks and labels from the main dotplot and overlay.
    # ax_overlay is independent (inset_axes) so set_xticks/yticks are safe here.
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax_overlay.set_xticks([])
    ax_overlay.set_yticks([])
    ax_overlay.xaxis.set_minor_locator(mticker.NullLocator())
    ax_overlay.yaxis.set_minor_locator(mticker.NullLocator())
    ax_overlay.set_xlabel("")
    ax_overlay.set_ylabel("")

    # Coordinate ticks on bar axes — set LAST because ax_xbar/ax_ybar share the major
    # tick container with ax/ax_top (sharex) and ax/ax_right (sharey).  Any earlier
    # set_xticks([]) call would wipe these, so they must come after all clearing.
    x_major, x_minor = _coord_ticks(asm2_segs_info, plot_w)
    ax_xbar.set_xticks(x_major)
    ax_xbar.xaxis.set_minor_locator(mticker.FixedLocator(x_minor))
    ax_xbar.tick_params(axis="x", which="major", bottom=True, top=False, labelbottom=False, length=5)
    ax_xbar.tick_params(axis="x", which="minor", bottom=True, top=False, labelbottom=False, length=3)
    # Segment labels at center of each segment, shifted down a bit below the bar.
    _xtrans = offset_copy(blended_transform_factory(ax_xbar.transData, ax_xbar.transAxes),
                          fig=fig, x=0, y=-4, units='points')
    for _, gen_s, gen_e, off in asm2_segs_info:
        ax_xbar.text(off + (gen_e - gen_s) / 2, 0, f"{gen_s:,}:{gen_e:,}",
                     transform=_xtrans, ha='center', va='top',
                     fontsize=7, rotation=45, clip_on=False)
    # Suppress both major and minor on all other x-sharing axes (minor container is also shared).
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=False)
    if ax_top is not None:
        ax_top.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=False)

    y_major, y_minor = _coord_ticks(asm1_segs_info, plot_h)
    ax_ybar.set_yticks(y_major)
    ax_ybar.yaxis.set_minor_locator(mticker.FixedLocator(y_minor))
    ax_ybar.tick_params(axis="y", which="major", left=True, right=False, labelleft=False, length=5)
    ax_ybar.tick_params(axis="y", which="minor", left=True, right=False, labelleft=False, length=3)
    # Segment labels at center of each segment, shifted right toward the bar.
    _ytrans = offset_copy(blended_transform_factory(ax_ybar.transAxes, ax_ybar.transData),
                          fig=fig, x=-8, y=0, units='points')
    for _, gen_s, gen_e, off in asm1_segs_info:
        ax_ybar.text(0, off + (gen_e - gen_s) / 2, f"{gen_s:,}:{gen_e:,}",
                     transform=_ytrans, ha='right', va='center',
                     fontsize=7, clip_on=False)
    # Suppress both major and minor on all other y-sharing axes.
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
    if ax_right is not None:
        ax_right.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)

    if title:
        fig.suptitle(title, fontsize=10, y=1.01)
    fig.savefig(output_path, bbox_inches="tight")
    _t0 = _t("savefig", _t0)
    if show:
        plt.show()
    plt.close(fig)
