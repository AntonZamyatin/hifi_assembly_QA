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

import pysam
import matplotlib.pyplot as plt
import math
import numpy as np
import sqlite3
import matplotlib.ticker as mticker
from aln_utils import cigar_stats
from typing import List
from SegmentStorage import SegmentCollection, Segment


PID_TIE_EPS_DEFAULT = 0.001
SUPPORT_COLORS = {
    "ref_only": "#2166ac",   # blue
    "ref_best": "#92c5de",   # light blue
    "cand_only": "#b2182b",  # red
    "cand_best": "#f4a3c2",  # pink
    "equal_best": "#9cd9af", # light green
}

def alignment_plot(bamfile: str, segments: List[Segment],
                   pid_threshold: float = 0.99, binsize: int = 1000,
                   ax=None, read_buckets=None, read_bucket_colors=None, orientation: str = "x",
                   return_max=False, pid_tie_eps: float = 1e-12):
    """
    Generate a simple alignment plot for a given BAM file, contig, and region.

    Parameters:
    - bamfile: Path to the BAM file containing alignments.
    - segments: List of Segment to plot alignments for (same assembly/contig).
    """
    standalone=False
    bam = pysam.AlignmentFile(bamfile, "rb")
    if not segments:
        raise ValueError("segments must be a non-empty list")
    segments = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    bins = [0.0 for i in range(math.ceil(sum((s.end - s.start) for s in segments) / binsize))]
    bucket_bins = None
    name_to_bucket = None
    if read_buckets is not None:
        bucket_bins = [[0.0 for _ in range(len(bins))] for _ in range(len(read_buckets))]
        name_to_bucket = {}
        for i, bucket in enumerate(read_buckets):
            for name in bucket:
                name_to_bucket[name] = i

    # Collect valid alignments first; later keep only best-PID alignments per read.
    valid_alns = []
    offset = 0
    for seg in segments:
        for aln in bam.fetch(seg.contig, start=seg.start, end=seg.end):
            pid = cigar_stats(aln)[-1]
            if pid >= pid_threshold:
                startbin = math.floor((max(aln.reference_start, seg.start) - seg.start + offset) / binsize)
                endbin = math.ceil((min(aln.reference_end, seg.end) - seg.start + offset) / binsize)
                if startbin < endbin:
                    valid_alns.append((aln.query_name, float(pid), startbin, endbin))
        offset += (seg.end - seg.start)

    # Keep only best-PID alignments for each read; split weight across ties.
    best_pid_by_read = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if read_name not in best_pid_by_read or pid > best_pid_by_read[read_name]:
            best_pid_by_read[read_name] = pid

    best_aln_counts_by_read = {}
    for read_name, pid, startbin, endbin in valid_alns:
        if pid + pid_tie_eps >= best_pid_by_read.get(read_name, float("-inf")):
            best_aln_counts_by_read[read_name] = best_aln_counts_by_read.get(read_name, 0) + 1

    for read_name, pid, startbin, endbin in valid_alns:
        if pid + pid_tie_eps < best_pid_by_read.get(read_name, float("-inf")):
            continue
        n = best_aln_counts_by_read.get(read_name, 1)
        w = 1.0 / float(max(n, 1))
        if bucket_bins is None:
            for i in range(startbin, endbin):
                bins[i] += w
        else:
            bucket_idx = name_to_bucket.get(read_name, None)
            if bucket_idx is None:
                for i in range(startbin, endbin):
                    bins[i] += w
            else:
                for i in range(startbin, endbin):
                    bucket_bins[bucket_idx][i] += w
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 5))
        standalone = True
    x = np.arange(len(bins)) * binsize
    if bucket_bins is None:
        if orientation == "x":
            ax.stackplot(x, bins, color='blue', step='pre')
        elif orientation == "y":
            ax.fill_betweenx(x, 0, bins, step='pre', color='blue')
        else:
            raise ValueError("orientation must be 'x' or 'y'")
    else:
        if read_bucket_colors is None:
            colors = [f"C{i % 10}" for i in range(len(bucket_bins))] + ['#9e9e9e']
        else:
            if len(read_bucket_colors) != len(bucket_bins):
                raise ValueError("read_bucket_colors length must match read_buckets length")
            colors = list(read_bucket_colors) + ['#6e6e6e']
        series = bucket_bins + [bins]
        if orientation == "x":
            ax.stackplot(x, series, colors=colors, step='pre')
        elif orientation == "y":
            y0 = np.zeros(len(bins))
            for s, c in zip(series, colors):
                y1 = y0 + np.array(s)
                ax.fill_betweenx(x, y0, y1, step='pre', color=c)
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
    #ax.set_ylim(bottom=0)

    if standalone:
        plt.xlabel('Position in Bin')
        plt.ylabel('Number of Alignments')
        plt.title('Alignment Plot for Segments (' + str(len(segments)) + ')')
        plt.show()
    if return_max:
        if bucket_bins is None:
            return ax, max(bins) if bins else 0.0
        total = np.array(bins, dtype=float)
        for s in bucket_bins:
            total += np.array(s, dtype=float)
        return ax, float(np.max(total)) if len(total) > 0 else 0.0
    return ax


def dotplot_kmer(seq_a: str, seq_b: str, k: int = 15, binsize: int = 1000,
                 ax=None, max_cells: int = 10_000_000,
                 mask_ambiguous: bool = True, log_scale: bool = False):
    """
    Generate a Gepard-like k-mer dot plot using binned k-mer matches.

    Parameters:
    - seq_a: First sequence string.
    - seq_b: Second sequence string.
    - k: k-mer size.
    - binsize: Window size for binning positions.
    - max_cells: Max allowed dense matrix cells (bins_a * bins_b).
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

    def _revcomp(seq: str) -> str:
        comp = str.maketrans("ACGTN", "TGCAN")
        return seq.translate(comp)[::-1]

    def _valid_kmer(kmer: str) -> bool:
        if not mask_ambiguous:
            return True
        for ch in kmer:
            if ch not in "ACGT":
                return False
        return True

    # Index k-mers in seq_a by bin.
    kmer_bins = {}
    for i in range(0, len_a - k + 1):
        kmer = seq_a[i:i + k]
        if not _valid_kmer(kmer):
            continue
        bin_a = i // binsize
        kmer_bins.setdefault(kmer, []).append(bin_a)
        kmer_bins.setdefault(_revcomp(kmer), []).append(bin_a)

    def _accumulate(seq: str, mat: np.ndarray):
        for j in range(0, len(seq) - k + 1):
            kmer = seq[j:j + k]
            if not _valid_kmer(kmer):
                continue
            bin_b = j // binsize
            bins = kmer_bins.get(kmer)
            if not bins:
                continue
            for bin_a in bins:
                mat[bin_a, bin_b] += 1

    # Accumulate shared k-mers per bin pair.
    mat = np.zeros((bins_a, bins_b), dtype=np.uint32)
    _accumulate(seq_b, mat)

    standalone = False
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))
        standalone = True

    plot_mat = np.log1p(mat) if log_scale else mat
    vmax = float(np.max(plot_mat))
    ax.imshow(plot_mat, origin='lower', aspect='equal', cmap='Greys', vmin=0,
              vmax=vmax)
    ax.set_xlabel("Seq B bins")
    ax.set_ylabel("Seq A bins")
    ax.set_title(f"K-mer dotplot (k={k}, binsize={binsize})")
    return ax


def reads_alignment_dotplot(sc: SegmentCollection, read_names: set[str], pid_threshold: float,
                            plen_threshold: float, db_path: str, ax=None, line_alpha: float = 1,
                            line_width_px: float = 6.0, pid_tie_eps: float = PID_TIE_EPS_DEFAULT):
    """
    Plot read support as a dot-plot (lines) between candidate (X) and reference (Y).
    Alignments are fetched from the read-indexed alignment database.
    """
    standalone = False if ax is not None else True
    def _build_axis_map(sc: List[Segment], assembly: str):
        segs = []
        offset = 0
        for seg in sc:
            if seg.assembly != assembly:
                continue
            segs.append((seg.contig, seg.start, seg.end, offset))
            offset += (seg.end - seg.start)
        return segs, offset

    def _map_pos(segs, contig: str, start: int, end: int):
        for c, s, e, off in segs:
            if c == contig and start >= s and end <= e:
                return off + (start - s), off + (end - s)
        return None

    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))

    sc = sorted(sc, key=lambda s: (s.contig, s.start, s.end))
    ref_segs, ref_len = _build_axis_map(sc, "ref")
    cand_segs, cand_len = _build_axis_map(sc, "cand")

    if not ref_segs or not cand_segs:
        raise ValueError("SegmentCollection must contain both ref and cand segments")

    def _fetch_alignments_local(con, reads, assemblies):
        if not reads:
            return {}
        asm_placeholders = ",".join(["?"] * len(assemblies))
        base_sql = f"""
            SELECT assembly, read_name, contig, start, end, strand, pid, plen
            FROM aln
            WHERE assembly IN ({asm_placeholders})
              AND contig IS NOT NULL
              AND start IS NOT NULL
              AND end IS NOT NULL
              AND read_name IN ({{READS_IN}})
        """
        out = {asm: [] for asm in assemblies}
        chunk_size = 800
        for i in range(0, len(reads), chunk_size):
            chunk = reads[i:i + chunk_size]
            reads_in = ",".join(["?"] * len(chunk))
            sql = base_sql.replace("{READS_IN}", reads_in)
            params = list(assemblies) + chunk
            rows = con.execute(sql, params).fetchall()
            for asm, read_name, contig, start, end, strand, pid, plen in rows:
                assert(pid is not None)
                assert(plen is not None)
                out[asm].append(
                    (str(read_name), str(contig), int(start), int(end),
                     int(strand), float(pid), float(plen))
                )
        return out

    # Fetch alignments for read names (includes strand).
    reads = list(read_names)
    if not reads:
        return ax
    con = sqlite3.connect(db_path)
    try:
        aln_by_asm = _fetch_alignments_local(con, reads, ["cand", "ref"])
    finally:
        con.close()

    # Group by read with pid and plen filter.
    by_read = {}
    for asm, alns in aln_by_asm.items():
        for read_name, contig, start, end, strand, pid, plen in alns:
            if pid < pid_threshold or plen < plen_threshold:
                continue
            by_read.setdefault(read_name, {}).setdefault(asm, []).append(
                (contig, start, end, strand, pid, plen)
            )

    # Build one plotted line per read using best-PID alignment on each assembly.
    lines = []
    for read_name in reads:
        per_asm = by_read.get(read_name, {})
        ref_alns = per_asm.get("ref", [])
        cand_alns = per_asm.get("cand", [])
        if not ref_alns or not cand_alns:
            continue
        r_contig, r_start, r_end, r_strand, r_pid, r_plen = max(
            ref_alns, key=lambda x: (x[4], x[5], x[2] - x[1])
        )
        c_contig, c_start, c_end, c_strand, c_pid, c_plen = max(
            cand_alns, key=lambda x: (x[4], x[5], x[2] - x[1])
        )

        r_mapped = _map_pos(ref_segs, r_contig, r_start, r_end)
        c_mapped = _map_pos(cand_segs, c_contig, c_start, c_end)
        if r_mapped is None or c_mapped is None:
            continue

        r0, r1 = r_mapped
        c0, c1 = c_mapped
        if r_strand == 0:
            r0, r1 = r1, r0
        if c_strand == 0:
            c0, c1 = c1, c0
        lines.append((c0, c1, r0, r1, c_pid, r_pid))

    # Plot fixed-width lines with PID-based colors.
    if lines:
        dpi = ax.figure.dpi if ax.figure is not None else 72.0
        lw_points = line_width_px * 72.0 / dpi
        for c0, c1, r0, r1, c_pid, r_pid in lines:
            if c_pid - r_pid > pid_tie_eps:
                color = SUPPORT_COLORS["cand_best"]
            elif r_pid - c_pid > pid_tie_eps:
                color = SUPPORT_COLORS["ref_best"]
            else:
                color = SUPPORT_COLORS["equal_best"]
            ax.plot([c0, c1], [r0, r1], color=color, alpha=line_alpha, linewidth=lw_points)

    # Draw segment boundaries (hardcoded styling)
    boundary_color = "#b0b0b0"
    boundary_alpha = 1
    boundary_lw = 1
    for _, start, end, off in cand_segs[:-1]:
        ax.axvline(off + (end - start), color=boundary_color,
                   alpha=boundary_alpha, linewidth=boundary_lw)
    for _, start, end, off in ref_segs[:-1]:
        ax.axhline(off + (end - start), color=boundary_color,
                   alpha=boundary_alpha, linewidth=boundary_lw)

    # Segment captions
    x_centers = [off + (end - start) / 2 for _, start, end, off in cand_segs]
    x_labels = [f"{contig}:{start}-{end}" for contig, start, end, _ in cand_segs]
    y_centers = [off + (end - start) / 2 for _, start, end, off in ref_segs]
    y_labels = [f"{contig}:{start}-{end}" for contig, start, end, _ in ref_segs]
    ax.set_xticks(x_centers)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_yticks(y_centers)
    ax.set_yticklabels(y_labels, rotation=0, ha="right", va="center")

    # Coordinate ticks every 50,000 bp (minor ticks)
    x_minor = np.arange(0, cand_len + 1, 50_000)
    y_minor = np.arange(0, ref_len + 1, 50_000)
    ax.xaxis.set_minor_locator(mticker.FixedLocator(x_minor))
    ax.yaxis.set_minor_locator(mticker.FixedLocator(y_minor))
    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

    # Styling: captions without tick lines, coordinate ticks with lines
    ax.tick_params(axis="x", which="major", length=0, pad=20)
    ax.tick_params(axis="y", which="major", length=0, pad=6)
    ax.tick_params(axis="x", which="minor", length=3, width=2, labelsize=8, labelbottom=True)
    ax.tick_params(axis="y", which="minor", length=3, width=2, labelleft=False)

    # Axis unit at end of X axis
    ax.text(1.0, -0.08, "Mbp", transform=ax.transAxes,
            ha="right", va="top")

    ax.set_xlabel("Candidate (concatenated segments)")
    ax.set_ylabel("Reference (concatenated segments)")
    if standalone:
        ax.set_xlim(0, max(cand_len - 1, 0))
        ax.set_ylim(0, max(ref_len - 1, 0))
        ax.set_aspect('equal', adjustable='box')
    return ax


def report_plot(seq_ref: str, seq_cand: str, segments: SegmentCollection, read_names: set[str],
                db_path: str, pid_threshold: float, plen_threshold: float, k: int = 15, binsize: int = 1000,
                output_path: str = "report_plot.png", dpi: int = 150,
                cand_bam: str | None = None, ref_bam: str | None = None):
    """
    Create a combined report plot with k-mer dotplot as background and
    read alignment dotplot overlaid on top. Saves to output_path.
    """

    def _total_len(sc, assembly: str) -> int:
        total = 0
        for seg in sc:
            if seg.assembly != assembly:
                continue
            total += (seg.end - seg.start)
        return total

    sc_sorted = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    cand_len = _total_len(sc_sorted, "cand")
    ref_len = _total_len(sc_sorted, "ref")

    height = max(4.0, 16.0 * (ref_len / max(cand_len, 1)) + 4)
    fig, ax = plt.subplots(figsize=(16.0, height), dpi=dpi)


    # Background k-mer dotplot
    dotplot_kmer(seq_ref, seq_cand, k=k, binsize=binsize, ax=ax, log_scale=True)

    if ax.images:
        ax.images[-1].set_extent((0, cand_len, 0, ref_len))


    # Full-size inset for reads overlay
    ax_overlay = ax.inset_axes((0, 0, 1, 1))
    ax_overlay.set_facecolor("none")
    ax_overlay.set_frame_on(False)
    reads_alignment_dotplot(segments, read_names, pid_threshold, plen_threshold, db_path, ax=ax_overlay)

    # Keep overlays aligned with base axes
    ax.set_xlim(0, cand_len)
    ax.set_ylim(0, ref_len)
    ax_overlay.set_xlim(ax.get_xlim())
    ax_overlay.set_ylim(ax.get_ylim())

    # Marginal alignment plots
    if cand_bam is not None and ref_bam is not None:
        cand_segs = [s for s in sc_sorted if s.assembly == "cand"]
        ref_segs = [s for s in sc_sorted if s.assembly == "ref"]

        def _interval_in_segments(seg_list: List[Segment], contig: str, start: int, end: int) -> bool:
            for seg in seg_list:
                if seg.contig == contig and start >= seg.start and end <= seg.end:
                    return True
            return False

        segs_by_asm = {"cand": cand_segs, "ref": ref_segs}

        def _classify_reads_for_support_colors(db_file: str, reads_set: set[str],
                                               pid_thr: float, eps: float = 0.001):
            reads = list(reads_set)
            out = {
                "ref_only": set(),          # blue
                "cand_only": set(),         # red
                "ref_best": set(),          # light blue
                "cand_best": set(),         # pink
                "equal_best": set(),        # light gray
            }
            if not reads:
                return out
            con = sqlite3.connect(db_file)
            try:
                chunk_size = 800
                for i in range(0, len(reads), chunk_size):
                    chunk = reads[i:i + chunk_size]
                    reads_in = ",".join(["?"] * len(chunk))
                    rows = con.execute(
                        f"""
                        SELECT assembly, read_name, contig, start, end, pid
                        FROM aln
                        WHERE assembly IN ('cand', 'ref')
                          AND read_name IN ({reads_in})
                        """,
                        chunk,
                    ).fetchall()
                    per_read = {}
                    for asm, rn, contig, start, end, pid in rows:
                        asm = str(asm)
                        if pid is None or float(pid) < pid_thr:
                            continue
                        if not _interval_in_segments(segs_by_asm[asm], str(contig), int(start), int(end)):
                            continue
                        if rn not in per_read:
                            per_read[rn] = {}
                        info = per_read[rn].setdefault(asm, {"n_aln": 0, "max_pid": float("-inf")})
                        info["n_aln"] += 1
                        info["max_pid"] = max(info["max_pid"], float(pid))
                    for rn in chunk:
                        asm_vals = per_read.get(rn, {})
                        ref_info = asm_vals.get("ref", {"n_aln": 0, "max_pid": 0.0})
                        cand_info = asm_vals.get("cand", {"n_aln": 0, "max_pid": 0.0})
                        has_ref = ref_info["n_aln"] > 0
                        has_cand = cand_info["n_aln"] > 0
                        if has_ref and not has_cand:
                            out["ref_only"].add(rn)
                        elif has_cand and not has_ref:
                            out["cand_only"].add(rn)
                        elif has_ref and has_cand:
                            ref_pid = ref_info["max_pid"]
                            cand_pid = cand_info["max_pid"]
                            d = cand_pid - ref_pid
                            if abs(d) <= eps:
                                out["equal_best"].add(rn)
                            elif d > 0:
                                out["cand_best"].add(rn)
                            else:
                                out["ref_best"].add(rn)
            finally:
                con.close()
            return out

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        ax_top = divider.append_axes("top", size=2, pad=0, sharex=ax)
        ax_right = divider.append_axes("right", size=2, pad=0, sharey=ax)
        ax_top.tick_params(axis="x", labelbottom=False)
        ax_right.tick_params(axis="y", labelleft=False)

        read_classes = _classify_reads_for_support_colors(
            db_path, read_names, pid_threshold, eps=PID_TIE_EPS_DEFAULT
        )
        cand_buckets = [
            read_classes["cand_only"],  # red
            read_classes["cand_best"],  # pink
            read_classes["equal_best"], # light gray
        ]
        cand_colors = [
            SUPPORT_COLORS["cand_only"],
            SUPPORT_COLORS["cand_best"],
            SUPPORT_COLORS["equal_best"],
        ]
        ref_buckets = [
            read_classes["ref_only"],   # blue
            read_classes["ref_best"],   # light blue
            read_classes["equal_best"], # light gray
        ]
        ref_colors = [
            SUPPORT_COLORS["ref_only"],
            SUPPORT_COLORS["ref_best"],
            SUPPORT_COLORS["equal_best"],
        ]
        ax_top, maxx = alignment_plot(cand_bam, cand_segs, pid_threshold=pid_threshold, binsize=binsize,
                                      ax=ax_top, read_buckets=cand_buckets,
                                      read_bucket_colors=cand_colors,
                                      orientation="x", return_max=True)
        ax_right, maxy = alignment_plot(ref_bam, ref_segs, pid_threshold=pid_threshold, binsize=binsize,
                                        ax=ax_right, read_buckets=ref_buckets,
                                        read_bucket_colors=ref_colors,
                                        orientation="y", return_max=True)

        m = max(maxx, maxy)
        ax_top.set_ylim(0, m)
        ax_right.set_xlim(0, m)
        ax_top.set_xticks([])
        ax_top.set_xlabel("")
        ax_right.set_yticks([])
        ax_right.set_ylabel("")

    # Remove ticks and labels
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax_overlay.set_xlabel("")
    ax_overlay.set_ylabel("")

    plt.show()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
