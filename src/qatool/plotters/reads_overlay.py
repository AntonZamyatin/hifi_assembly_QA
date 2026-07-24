"""Read-support dot plot: one line per read between asm2 (X) and asm1 (Y)."""

import sqlite3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.collections import LineCollection

from qatool.SegmentStorage import SegmentCollection
from qatool.plotters.axis import _build_axis_map
from qatool.plotters.colors import SUPPORT_COLORS, PID_TIE_EPS_DEFAULT
from qatool.plotters.db import _fetch_aln_from_db


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
