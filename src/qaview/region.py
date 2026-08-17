"""RegionData — everything one interactive view needs, as plain arrays.

No matplotlib artists and no Bokeh objects: this layer is renderer-agnostic and
JSON-serialisable, so the same numbers feed the interactive page and the published
`report_plot` figure.

Coordinates here are **exact bp** (`_build_axis_map` with binsize=None), unlike the
bin-padded space `report_plot` uses.  The padding exists only so k-mer mosaic pixels
never straddle a segment boundary, and the interactive view draws no mosaic by
default.  `ViewState` will therefore carry genomic coordinates, not plot ones, so
each renderer can build its own mapping — see `docs/qaview.md` §5.
"""

import math
import sqlite3
from dataclasses import dataclass

import numpy as np
import pysam

from qatool.plotters.axis import _build_axis_map, plot_to_genome
from qatool.plotters.classify import _classify_reads_for_support_colors
from qatool.plotters.colors import BEST_PID_DELTA_DEFAULT, PID_TIE_EPS_DEFAULT
from qatool.plotters.coverage import CoverageBins, compute_coverage_bins
from qatool.plotters.db import _collect_bam_read_names, _fetch_aln_from_db
from qatool.plotters.reads_dotplot import (
    _build_lookup, _fetch_aln_blocks, _map_in_segment, _read_span, _read_to_ref,
    _visible_read_span,
)

# Bucket order used for the stacked marginal coverage plots; mirrors report.py so the
# two renderers stack the same categories in the same order.
SUPPORT_BUCKET_KEYS = ["asm1_only", "asm2_only", "asm1_best", "asm2_best", "common", "other"]

# Painter's-algorithm layers for the dotplot.  Rows are stored in this order so a
# single glyph renderer draws them back-to-front — which keeps line row indices
# consistent with the CSR highlight indexes, unlike splitting into four renderers.
LAYER_COMMON, LAYER_OTHER, LAYER_SUPPORT, LAYER_SUPPLEMENTARY = 0, 1, 2, 3

_SUPPORT_CLASSES = {"asm1_only", "asm2_only", "asm1_best", "asm2_best"}


def _layer_of(read_class: str, supplementary: bool) -> int:
    """Which z-layer a line belongs to.  Supplementary wins over classification."""
    if supplementary:
        return LAYER_SUPPLEMENTARY
    if read_class in _SUPPORT_CLASSES:
        return LAYER_SUPPORT
    if read_class == "common":
        return LAYER_COMMON
    return LAYER_OTHER


@dataclass
class ReadLines:
    """Read-correspondence line segments in plot space (asm2 = X, asm1 = Y).

    One entry per (asm1 block, asm2 block) pair sharing read span.  Unlike
    `reads_block_dotplot` these are **not** merged: merging is resolution-dependent
    and would have to be redone on every zoom, so the interactive renderer encodes
    depth with alpha instead.

    The `contig*` / `ref*` / `pid*` columns describe the **drawn** (clipped) extent,
    not the full alignment, so a tooltip reports what is actually on screen.
    """
    x0: np.ndarray
    y0: np.ndarray
    x1: np.ndarray
    y1: np.ndarray
    read_idx: np.ndarray
    supplementary: np.ndarray
    n_dropped: int
    layer: np.ndarray
    contig1: list
    contig2: list
    ref1_0: np.ndarray
    ref1_1: np.ndarray
    ref2_0: np.ndarray
    ref2_1: np.ndarray
    pid1: np.ndarray
    pid2: np.ndarray

    def __len__(self) -> int:
        return len(self.x0)

    def by_read(self, n_reads: int) -> tuple:
        """CSR map read index -> line row indices, for client-side highlighting.

        Returns (offsets, indices) with offsets of length n_reads + 1.
        """
        counts = np.zeros(n_reads, dtype=np.int64)
        for ri in self.read_idx:
            counts[ri] += 1
        offsets = np.zeros(n_reads + 1, dtype=np.int64)
        np.cumsum(counts, out=offsets[1:])
        cursor = offsets[:-1].copy()
        indices = np.zeros(len(self.read_idx), dtype=np.int32)
        for row, ri in enumerate(self.read_idx):
            indices[cursor[ri]] = row
            cursor[ri] += 1
        return offsets, indices


@dataclass
class RegionData:
    """One recruited region, reduced to arrays."""
    asm1: str
    asm2: str
    binsize: int
    segments: list
    asm1_segs_info: list
    asm2_segs_info: list
    plot_w: int                 # asm2 extent (X), exact bp
    plot_h: int                 # asm1 extent (Y), exact bp
    read_names: list
    read_class: dict
    lines: ReadLines
    cov1: CoverageBins          # asm1 — right-hand marginal
    cov2: CoverageBins          # asm2 — top marginal
    bucket_keys: list

    def segs_info(self, assembly: str) -> list:
        if assembly == self.asm1:
            return self.asm1_segs_info
        if assembly == self.asm2:
            return self.asm2_segs_info
        raise ValueError(f"unknown assembly {assembly!r}")

    def genome_at(self, assembly: str, pos: float):
        """Plot coordinate -> (contig, genomic position) for the crosshair readout."""
        return plot_to_genome(self.segs_info(assembly), pos)

    def bin_edges(self, assembly: str) -> tuple:
        """Exact-bp (left, right) edges for each coverage bin of `assembly`.

        Matches the bin order of the corresponding CoverageBins: segments in
        (contig, start, end) order, ceil(seg_len / binsize) bins each, final bin of a
        segment clipped to the segment end.
        """
        left, right = [], []
        for _, start, end, offset in self.segs_info(assembly):
            seg_len = end - start
            for j in range(math.ceil(seg_len / self.binsize)):
                l = offset + j * self.binsize
                left.append(l)
                right.append(min(offset + seg_len, l + self.binsize))
        return np.array(left, dtype=np.int64), np.array(right, dtype=np.int64)

    @classmethod
    def build(cls, segments, read_names, db_path: str,
              asm1_bam: str, asm2_bam: str,
              asm1_fasta: str | None = None, asm2_fasta: str | None = None,
              asm1: str = "asm1", asm2: str = "asm2",
              pid_threshold: float = 0.99, plen_threshold: float = 0.99,
              binsize: int | None = None,
              best_pid_delta: float = BEST_PID_DELTA_DEFAULT,
              pid_tie_eps: float = PID_TIE_EPS_DEFAULT) -> "RegionData":
        """Assemble a RegionData from a recruited SegmentCollection and read set.

        binsize=None auto-sizes to ~1000 bins on the longer axis, the same rule
        `report_plot._resolve_binsize` uses (minus its k-mer floor, which does not
        apply here — no mosaic is drawn).
        """
        sc_sorted = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
        if binsize is None:
            total_a = sum(s.end - s.start for s in sc_sorted if s.assembly == asm1)
            total_b = sum(s.end - s.start for s in sc_sorted if s.assembly == asm2)
            binsize = max(1, math.ceil(max(total_a, total_b, 1) / 1000))
        # Sorting once here keeps _build_axis_map's iteration order and
        # compute_coverage_bins' internal sort in agreement, so bin i of the coverage
        # arrays and bin i of bin_edges() refer to the same piece of genome.
        asm1_segs_info, plot_h = _build_axis_map(sc_sorted, asm1)
        asm2_segs_info, plot_w = _build_axis_map(sc_sorted, asm2)
        asm1_segs = [s for s in sc_sorted if s.assembly == asm1]
        asm2_segs = [s for s in sc_sorted if s.assembly == asm2]
        if not asm1_segs or not asm2_segs:
            raise ValueError("segments must contain both asm1 and asm2 segments")

        # Expand the seed set so classification covers every read visible in a marginal.
        all_reads = set(read_names)
        all_reads |= _collect_bam_read_names(asm2_bam, asm2_segs, pid_threshold, plen_threshold)
        all_reads |= _collect_bam_read_names(asm1_bam, asm1_segs, pid_threshold, plen_threshold)

        con = sqlite3.connect(db_path)
        try:
            aln_by_asm = _fetch_aln_from_db(con, list(all_reads), [asm1, asm2],
                                            pid_threshold, plen_threshold)
            aln_blocks = _fetch_aln_blocks(con, list(all_reads), [asm1, asm2],
                                           pid_threshold, plen_threshold)
        finally:
            con.close()

        read_classes = _classify_support(db_path, all_reads, pid_threshold, plen_threshold,
                                         asm1, asm2, pid_tie_eps, aln_by_asm,
                                         asm1_fasta, asm2_fasta, best_pid_delta)
        read_names_l = sorted(all_reads)
        read_pos = {n: i for i, n in enumerate(read_names_l)}
        read_class = {n: k for k in SUPPORT_BUCKET_KEYS for n in read_classes[k]}

        lines = _compute_lines(aln_blocks, read_names_l, read_pos,
                               asm1_segs_info, asm2_segs_info, asm1, asm2)
        lines = _order_for_drawing(lines, read_names_l, read_class)

        buckets = [read_classes[k] for k in SUPPORT_BUCKET_KEYS]
        cov1 = compute_coverage_bins(asm1_bam, asm1_segs, pid_threshold, plen_threshold,
                                     binsize, buckets, best_pid_delta)
        cov2 = compute_coverage_bins(asm2_bam, asm2_segs, pid_threshold, plen_threshold,
                                     binsize, buckets, best_pid_delta)

        return cls(asm1=asm1, asm2=asm2, binsize=binsize, segments=sc_sorted,
                   asm1_segs_info=asm1_segs_info, asm2_segs_info=asm2_segs_info,
                   plot_w=plot_w, plot_h=plot_h,
                   read_names=read_names_l, read_class=read_class, lines=lines,
                   cov1=cov1, cov2=cov2, bucket_keys=list(SUPPORT_BUCKET_KEYS))


def _classify_support(db_path, all_reads, pid_threshold, plen_threshold,
                      asm1, asm2, pid_tie_eps, aln_by_asm,
                      asm1_fasta, asm2_fasta, best_pid_delta) -> dict:
    """Support-classify every read globally (segs_by_asm=None), as report_plot does,
    so a zoomed subset keeps each read's globally-correct category."""
    fa1 = pysam.FastaFile(asm1_fasta) if asm1_fasta else None
    fa2 = pysam.FastaFile(asm2_fasta) if asm2_fasta else None
    try:
        return _classify_reads_for_support_colors(
            db_path, all_reads, pid_threshold, plen_threshold,
            segs_by_asm=None, asm1=asm1, asm2=asm2, pid_tie_eps=pid_tie_eps,
            aln_by_asm=aln_by_asm, asm1_fasta=fa1, asm2_fasta=fa2,
            best_pid_delta=best_pid_delta,
        )
    finally:
        if fa1:
            fa1.close()
        if fa2:
            fa2.close()


def _compute_lines(aln_blocks, read_names, read_pos,
                   asm1_segs_info, asm2_segs_info, asm1, asm2) -> ReadLines:
    """Thread each read through its asm1 and asm2 blocks into plot-space line segments.

    Geometry helpers are reused verbatim from `qatool.plotters.reads_dotplot`; only
    the merging and matplotlib drawing are left behind.
    """
    y_lookup = _build_lookup(asm1_segs_info)   # asm1 -> Y
    x_lookup = _build_lookup(asm2_segs_info)   # asm2 -> X

    by_read: dict = {}
    for asm, blocks in aln_blocks.items():
        for b in blocks:
            by_read.setdefault(b[0], {}).setdefault(asm, []).append(b)

    x0s, y0s, x1s, y1s, idxs, supps = [], [], [], [], [], []
    ct1s, ct2s, r1_0, r1_1, r2_0, r2_1, p1s, p2s = [], [], [], [], [], [], [], []
    n_dropped = 0
    for rn in read_names:
        per = by_read.get(rn, {})
        a1 = per.get(asm1, [])
        a2 = per.get(asm2, [])
        if not a1 or not a2:
            continue
        ri = read_pos[rn]
        for (_, c1, s1, e1, st1, sl1, sr1, hl1, hr1, q1, t1, pid1) in a1:
            r0a, r1a = _read_span(st1, sl1, sr1, hl1, hr1, q1)
            v1 = _visible_read_span(y_lookup, c1, s1, e1, r0a, r1a, st1)
            if v1 is None:
                n_dropped += 1
                continue                  # block lands on no displayed segment at all
            for (_, c2, s2, e2, st2, sl2, sr2, hl2, hr2, q2, t2, pid2) in a2:
                r0b, r1b = _read_span(st2, sl2, sr2, hl2, hr2, q2)
                v2 = _visible_read_span(x_lookup, c2, s2, e2, r0b, r1b, st2)
                if v2 is None:
                    n_dropped += 1
                    continue
                # Blocks overhanging a segment edge are clipped to their visible part
                # rather than discarded; n_dropped now counts only genuine misses.
                ro0 = max(r0a, r0b, v1[0], v2[0])
                ro1 = min(r1a, r1b, v1[1], v2[1])
                if ro1 <= ro0:
                    n_dropped += 1
                    continue              # no shared, visible stretch of read
                ra0 = _read_to_ref(ro0, r0a, r1a, s1, e1, st1)
                ra1 = _read_to_ref(ro1, r0a, r1a, s1, e1, st1)
                rb0 = _read_to_ref(ro0, r0b, r1b, s2, e2, st2)
                rb1 = _read_to_ref(ro1, r0b, r1b, s2, e2, st2)
                x0s.append(_map_in_segment(v2[2], rb0))
                y0s.append(_map_in_segment(v1[2], ra0))
                x1s.append(_map_in_segment(v2[2], rb1))
                y1s.append(_map_in_segment(v1[2], ra1))
                idxs.append(ri)
                supps.append((t1 == "supplementary") or (t2 == "supplementary"))
                ct1s.append(c1); ct2s.append(c2)
                r1_0.append(min(ra0, ra1)); r1_1.append(max(ra0, ra1))
                r2_0.append(min(rb0, rb1)); r2_1.append(max(rb0, rb1))
                p1s.append(pid1); p2s.append(pid2)

    return ReadLines(
        x0=np.array(x0s, dtype=float), y0=np.array(y0s, dtype=float),
        x1=np.array(x1s, dtype=float), y1=np.array(y1s, dtype=float),
        read_idx=np.array(idxs, dtype=np.int32),
        supplementary=np.array(supps, dtype=bool),
        n_dropped=n_dropped,
        layer=np.zeros(len(x0s), dtype=np.int8),
        contig1=ct1s, contig2=ct2s,
        ref1_0=np.array(r1_0, dtype=np.int64), ref1_1=np.array(r1_1, dtype=np.int64),
        ref2_0=np.array(r2_0, dtype=np.int64), ref2_1=np.array(r2_1, dtype=np.int64),
        pid1=np.array(p1s, dtype=float), pid2=np.array(p2s, dtype=float),
    )


def _order_for_drawing(lines: ReadLines, read_names: list, read_class: dict) -> ReadLines:
    """Reorder rows back-to-front: common, unclassified, support, supplementary.

    Sorting the stored rows (rather than splitting into per-layer renderers) means the
    dotplot stays one glyph with one data source, so `by_read()` row indices remain
    valid for the hover/selection callbacks.  `kind="stable"` keeps rows of the same
    layer in their original read order.
    """
    if len(lines) == 0:
        return lines
    layer = np.array([
        _layer_of(read_class.get(read_names[ri], "other"), bool(sup))
        for ri, sup in zip(lines.read_idx, lines.supplementary)
    ], dtype=np.int8)
    order = np.argsort(layer, kind="stable")
    return ReadLines(
        x0=lines.x0[order], y0=lines.y0[order], x1=lines.x1[order], y1=lines.y1[order],
        read_idx=lines.read_idx[order], supplementary=lines.supplementary[order],
        n_dropped=lines.n_dropped, layer=layer[order],
        contig1=[lines.contig1[i] for i in order], contig2=[lines.contig2[i] for i in order],
        ref1_0=lines.ref1_0[order], ref1_1=lines.ref1_1[order],
        ref2_0=lines.ref2_0[order], ref2_1=lines.ref2_1[order],
        pid1=lines.pid1[order], pid2=lines.pid2[order],
    )
