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
from dataclasses import dataclass, field, replace

import numpy as np
import pysam

from qatool.plotters.axis import _build_axis_map
from qatool.plotters.classify import _classify_reads_for_support_colors
from qatool.plotters.colors import BEST_PID_DELTA_DEFAULT, PID_TIE_EPS_DEFAULT
from qatool.plotters.coverage import CoverageBins
from qatool.plotters.db import _collect_bam_read_names, _fetch_aln_from_db
from qatool.plotters.reads_dotplot import (
    _BLOCK_COLS, _build_lookup, _fetch_aln_blocks, _map_in_segment, _map_point, _read_span,
    _read_to_ref, _visible_read_span,
)

# Bucket order used for the stacked marginal coverage plots; mirrors report.py so the
# two renderers stack the same categories in the same order.
SUPPORT_BUCKET_KEYS = ["asm1_only", "asm2_only", "asm1_best", "asm2_best", "common", "other"]

# Painter's-algorithm layers for the dotplot.  Rows are stored in this order so a
# single glyph renderer draws them back-to-front — which keeps line row indices
# consistent with the CSR highlight indexes, unlike splitting into four renderers.
LAYER_COMMON, LAYER_OTHER, LAYER_SUPPORT, LAYER_SUPPLEMENTARY = 0, 1, 2, 3
# A one-sided read's partial alignments in the other candidate: never drawn but for
# the reads a hover or a tap selects (`_partial_lines`).
LAYER_PARTIAL = 4

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
    # per row, what a partial line is ("partial supplementary in curated, ..."); "" for
    # an ordinary row; empty when there are no partial rows
    part: list = field(default_factory=list)

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
class Placements:
    """Each read's best placements in one assembly, as intervals on its axis.

    One row per (placement, displayed segment it overlaps): [p0, p1) in plot
    coordinates, weight 1/n (n placements genome-wide), the read's index into
    `RegionData.read_names`, and the segment's index into the axis map. Rows are
    sorted by read, so read g's rows are `first[g]:last[g]`. The marginals bin these
    afresh for whatever range is in view (`marginal_bins`, glyphs' REBIN_JS).
    """
    p0: np.ndarray
    p1: np.ndarray
    w: np.ndarray
    read: np.ndarray
    seg: np.ndarray
    first: np.ndarray
    last: np.ndarray

    def __len__(self) -> int:
        return len(self.p0)


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
    # assembly -> {contig: length} for the contigs in `segments`, from the BAM
    # headers: which segment edges are contig ends (glyphs.contig_bar_figure).
    contig_lengths: dict = field(default_factory=dict)
    # rows of `lines` that are partial alignments (LAYER_PARTIAL), for counting
    n_partial: int = 0
    # large indels on those rows, shown with them: x, y, read (index), text
    indels: dict = field(default_factory=dict)
    # assembly -> Placements: what the marginals are binned from, at any bin size
    placements: dict = field(default_factory=dict)
    # assembly -> one flag per segs_info row: the segment runs backwards along the
    # axis (`arrange`); missing means none does
    flips: dict = field(default_factory=dict)

    def segs_info(self, assembly: str) -> list:
        if assembly == self.asm1:
            return self.asm1_segs_info
        if assembly == self.asm2:
            return self.asm2_segs_info
        raise ValueError(f"unknown assembly {assembly!r}")

    def flips_of(self, assembly: str) -> list:
        """Per segs_info row: is that segment drawn backwards?"""
        return list(self.flips.get(assembly) or [False] * len(self.segs_info(assembly)))

    def genome_at(self, assembly: str, pos: float):
        """Plot coordinate -> (contig, genomic position) for the crosshair readout."""
        for (contig, start, end, off), f in zip(self.segs_info(assembly), self.flips_of(assembly)):
            if off <= pos < off + (end - start):
                return contig, (end - 1 - int(pos - off)) if f else start + int(pos - off)
        return None

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
              pid_tie_eps: float = PID_TIE_EPS_DEFAULT,
              supp_plen_threshold: float | None = None,
              partial_pid_threshold: float | None = None,
              min_indel: int = 50) -> "RegionData":
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

        # Every read is shown at its best placements only, whatever its class:
        # in the dotplot it runs between its best placements in the two
        # assemblies, and in each marginal its weight of 1 is split evenly over
        # its best placements in that assembly, counted over the whole assembly
        # (`_best_placements`). The aln rows are fetched by read name, genome-wide.
        best = {asm: _best_by_read(aln_blocks[asm], best_pid_delta) for asm in (asm1, asm2)}
        lines = _compute_lines(best, read_names_l, read_pos,
                               asm1_segs_info, asm2_segs_info, asm1, asm2)
        lines = _order_for_drawing(lines, read_names_l, read_class)

        # A read only one candidate holds whole may be a breakpoint read: its
        # partial alignments in the other candidate inside the locus are added as
        # hidden rows, drawn when a hover or a tap selects the read. They change
        # nothing else -- not its class, the marginals, or the default view.
        n_partial, indels = 0, {}
        if supp_plen_threshold is not None:
            only = {asm1: read_classes["asm2_only"], asm2: read_classes["asm1_only"]}
            lengths = {asm1: _all_contig_lengths(asm1_bam), asm2: _all_contig_lengths(asm2_bam)}
            ppid = pid_threshold if partial_pid_threshold is None else partial_pid_threshold
            pieces = _partial_alignments(db_path, only, pid_threshold, plen_threshold,
                                         supp_plen_threshold, lengths, ppid, min_indel,
                                         best_pid_delta)
            plines, indels = _partial_lines(pieces, best, read_pos, asm1_segs_info,
                                            asm2_segs_info, asm1, asm2)
            n_partial = len(plines)
            lines = _concat_lines(lines, plines)

        buckets = [read_classes[k] for k in SUPPORT_BUCKET_KEYS]
        cov1 = _coverage_at_best(best[asm1], asm1_segs, binsize, buckets)
        cov2 = _coverage_at_best(best[asm2], asm2_segs, binsize, buckets)

        return cls(asm1=asm1, asm2=asm2, binsize=binsize, segments=sc_sorted,
                   asm1_segs_info=asm1_segs_info, asm2_segs_info=asm2_segs_info,
                   plot_w=plot_w, plot_h=plot_h,
                   read_names=read_names_l, read_class=read_class, lines=lines,
                   cov1=cov1, cov2=cov2, bucket_keys=list(SUPPORT_BUCKET_KEYS),
                   contig_lengths={asm1: _contig_lengths(asm1_bam, asm1_segs),
                                   asm2: _contig_lengths(asm2_bam, asm2_segs)},
                   n_partial=n_partial, indels=indels,
                   placements={asm1: _placements(best[asm1], asm1_segs_info, read_pos),
                               asm2: _placements(best[asm2], asm2_segs_info, read_pos)})


def _contig_lengths(bam_path: str, segs) -> dict:
    """{contig: length} for the segments' contigs, from the BAM header."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        known = set(bam.references)
        return {c: bam.get_reference_length(c) for c in {s.contig for s in segs} if c in known}


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


_PRIMARY = ("primary", "primary_wsup")


# A partial alignment whose clipped side ends this close to its contig's end runs
# off the contig: a contig break rather than a breakpoint inside it.
CLIP_END_TOL = 1000
_PARTIAL_KIND = {"primary": "primary, clipped", "primary_wsup": "primary, split",
                 "supplementary": "supplementary", "secondary": "secondary"}


def _all_contig_lengths(bam_path: str) -> dict:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return dict(zip(bam.references, bam.lengths))


_CIGAR = __import__("re").compile(r"(\d+)([MIDNSHP=X])")
# At most this many of a read's other alignments are kept: in a repeat a read can
# have dozens tied at its best, each drawn against each of its whole placements.
PARTIAL_MAX_PER_READ = 5


def _large_indels(cigar: str, start: int, min_len: int) -> list:
    """(reference position, "insertion" | "deletion", length) for each indel of at
    least `min_len` bp in an alignment starting at `start`; a deletion is placed at
    the middle of the reference it skips, an insertion where it occurs."""
    out, ref = [], start
    for n, op in _CIGAR.findall(cigar or ""):
        n = int(n)
        if op == "D" and n >= min_len:
            out.append((ref + n // 2, "deletion", n))
        elif op == "I" and n >= min_len:
            out.append((ref, "insertion", n))
        if op in "MD=XN":
            ref += n
    return out


def _partial_alignments(db_path: str, reads_by_asm: dict, pid: float, plen: float,
                        min_plen: float, lengths: dict, partial_pid: float | None = None,
                        min_indel: int = 50, best_pid_delta: float = BEST_PID_DELTA_DEFAULT,
                        chunk: int = 800) -> dict:
    """assembly -> {read: [(row, label, indels)]}: the reads' alignments there that do
    not qualify at the display pid and plen, but reach `partial_pid` (default the
    display pid) with plen at least `min_plen` -- partial ones of any type, and
    whole-length ones short of the display pid that hold an indel of `min_indel` bp
    or more: an SV kept inside one alignment as a single long gap (one without is
    divergence, not an SV, and is left out). `indels` are those of at least `min_indel` bp (`_large_indels`).
    Only a read's best of these are kept, as for placements: non-supplementary ones
    within `best_pid_delta` of its best pid among them (at most
    PARTIAL_MAX_PER_READ), and the supplementary pieces of a kept primary -- a split
    read keeps its pieces, a read in a repeat not every copy.
    Rows as _fetch_aln_blocks (sl, sr, hl, hr at 5-8)."""
    ppid = pid if partial_pid is None else partial_pid
    out: dict = {}
    con = sqlite3.connect(db_path)
    try:
        for asm, reads in reads_by_asm.items():
            reads = list(reads)
            got: dict = {}
            for i in range(0, len(reads), chunk):
                c = reads[i:i + chunk]
                rows = con.execute(
                    f"SELECT {_BLOCK_COLS}, plen, cigar FROM aln WHERE assembly = ? "
                    f"AND contig IS NOT NULL AND pid >= ? AND plen >= ? "
                    f"AND NOT (pid >= ? AND plen >= ?) "
                    f"AND read_name IN ({','.join('?' * len(c))})",
                    [asm, ppid, min_plen, pid, plen, *c])
                for r in rows:
                    n = lengths[asm].get(r[1], 0)
                    at_end = ((r[5] + r[7] > 0 and r[2] <= CLIP_END_TOL)
                              or (r[6] + r[8] > 0 and r[3] >= n - CLIP_END_TOL))
                    whole = r[12] >= plen
                    big = _large_indels(r[13], r[2], min_indel)
                    if whole and not big:
                        continue      # short of the pid throughout: divergence, not an SV
                    kind = "whole-length" if whole else "partial " + _PARTIAL_KIND.get(r[10], r[10])
                    label = (f"{kind} in {asm}, pid {r[11]:.4f}, plen {r[12]:.2f}"
                             + ("" if whole else f", {'at a contig end' if at_end else 'inside a contig'}"))
                    if big:
                        label += "; " + ", ".join(f"{ln:,} bp {k} at {r[1]}:{pos + 1:,}"
                                                  for pos, k, ln in big[:3])
                        if len(big) > 3:
                            label += f" (+{len(big) - 3} more)"
                    got.setdefault(r[0], []).append((r[:12], label, big))
            out[asm] = {rn: _best_partials(items, best_pid_delta) for rn, items in got.items()}
    finally:
        con.close()
    return out


def _best_partials(items: list, tol: float) -> list:
    """A read's best other alignments (see `_partial_alignments`): rows are
    (row, label, indels), row[10] the type and row[11] the pid."""
    main = sorted((x for x in items if x[0][10] != "supplementary"), key=lambda x: -x[0][11])
    if not main:
        return []
    kept = [x for x in main if x[0][11] + tol >= main[0][0][11]][:PARTIAL_MAX_PER_READ]
    if any(x[0][10] in _PRIMARY for x in kept):
        kept += [x for x in items if x[0][10] == "supplementary"]
    return kept


def _partial_lines(pieces: dict, best: dict, read_pos: dict, asm1_segs_info, asm2_segs_info,
                   asm1: str, asm2: str) -> tuple:
    """(lines, indels): lines from each one-sided read's best placements in the
    candidate holding it whole to each of its other alignments in the other one,
    LAYER_PARTIAL, labelled per row; and a point on those lines for each large indel
    -- x, y, read (index into the region's reads), text."""
    parts = []
    marks = {"x": [], "y": [], "read": [], "text": []}
    lookup = {asm1: _build_lookup(asm1_segs_info), asm2: _build_lookup(asm2_segs_info)}
    for asm, by_read in pieces.items():
        whole = asm2 if asm == asm1 else asm1
        for rn, items in sorted(by_read.items()):
            full = best[whole].get(rn)
            if not full or not full[1]:
                continue
            for row, label, big in items:
                side = {whole: {rn: full}, asm: {rn: ([row], 1)}}
                ln = _compute_lines(side, [rn], read_pos, asm1_segs_info, asm2_segs_info, asm1, asm2)
                if not len(ln):
                    continue
                ln.layer = np.full(len(ln), LAYER_PARTIAL, dtype=np.int8)
                ln.part = [label] * len(ln)
                parts.append(ln)
                for pos, k, n in big:
                    pt = _map_point(lookup[asm], row[1], pos, pos, pos)
                    if pt is None:
                        continue
                    for i in range(len(ln)):
                        a0, a1 = (ln.x0[i], ln.x1[i]) if asm == asm2 else (ln.y0[i], ln.y1[i])
                        if a0 == a1 or not min(a0, a1) <= pt <= max(a0, a1):
                            continue
                        t = (pt - a0) / (a1 - a0)
                        marks["x"].append(float(ln.x0[i] + t * (ln.x1[i] - ln.x0[i])))
                        marks["y"].append(float(ln.y0[i] + t * (ln.y1[i] - ln.y0[i])))
                        marks["read"].append(read_pos[rn])
                        marks["text"].append(f"{n:,} bp {k} in {asm} at {row[1]}:{pos + 1:,}")
    out = _compute_lines({asm1: {}, asm2: {}}, [], read_pos, asm1_segs_info, asm2_segs_info,
                         asm1, asm2)
    for ln in parts:
        out = _concat_lines(out, ln)
    return out, (marks if marks["x"] else {})


def _concat_lines(a: "ReadLines", b: "ReadLines") -> "ReadLines":
    """Rows of `a` then `b`; `part` filled with "" where a side has none."""
    if not len(b):
        return a
    cat = lambda x, y: np.concatenate([x, y])
    return ReadLines(
        x0=cat(a.x0, b.x0), y0=cat(a.y0, b.y0), x1=cat(a.x1, b.x1), y1=cat(a.y1, b.y1),
        read_idx=cat(a.read_idx, b.read_idx).astype(np.int32),
        supplementary=cat(a.supplementary, b.supplementary).astype(bool),
        n_dropped=a.n_dropped + b.n_dropped, layer=cat(a.layer, b.layer).astype(np.int8),
        contig1=list(a.contig1) + list(b.contig1), contig2=list(a.contig2) + list(b.contig2),
        ref1_0=cat(a.ref1_0, b.ref1_0), ref1_1=cat(a.ref1_1, b.ref1_1),
        ref2_0=cat(a.ref2_0, b.ref2_0), ref2_1=cat(a.ref2_1, b.ref2_1),
        pid1=cat(a.pid1, b.pid1), pid2=cat(a.pid2, b.pid2),
        part=(list(a.part) or [""] * len(a)) + (list(b.part) or [""] * len(b)),
    )


def _best_placements(blocks: list, tol: float) -> tuple:
    """(kept, n): a read's alignments at its best placements in one assembly.

    Its non-supplementary alignments within `tol` (best_pid_delta) of its best pid
    there, plus the supplementary pieces of a best primary (they belong to it, at
    lower pid); `n` counts the non-supplementary ones -- the placements its weight
    of 1 is split over. `blocks` are all its qualifying alignments in the assembly,
    genome-wide; a read with only supplementary ones there has none.
    """
    main = [b for b in blocks if b[10] != "supplementary"]
    if not main:
        return [], 0
    best = max(b[11] for b in main)
    kept = [b for b in main if b[11] + tol >= best]
    if any(b[10] in _PRIMARY for b in kept):
        kept += [b for b in blocks if b[10] == "supplementary"]
    return kept, sum(b[10] != "supplementary" for b in kept)


def read_lengths(db_path: str, names: list, assemblies: list, chunk: int = 800) -> dict:
    """read -> its length (q_len) from the aln table."""
    out = {}
    ph = ",".join("?" * len(assemblies))
    con = sqlite3.connect(db_path)
    try:
        for i in range(0, len(names), chunk):
            c = list(names[i:i + chunk])
            sql = (f"SELECT read_name, MAX(q_len) FROM aln WHERE assembly IN ({ph}) "
                   f"AND read_name IN ({','.join('?' * len(c))}) GROUP BY read_name")
            out.update(con.execute(sql, list(assemblies) + c).fetchall())
    finally:
        con.close()
    return out


def read_alignments(db_path: str, name: str, assemblies: list, pid: float, plen: float,
                    tol: float) -> list:
    """Every alignment of read `name` in `assemblies`, genome-wide, as dicts:
    assembly, contig, start, end, strand ("+"/"-"), pid, plen, type, and two marks --
    `passes` (the display pid/plen) and `best` (one of its best placements there,
    `_best_placements`, the ones the plots draw). Unmapped records are left out."""
    con = sqlite3.connect(db_path)
    try:
        rows = con.execute(
            f"SELECT assembly, {_BLOCK_COLS}, plen FROM aln WHERE assembly IN "
            f"({','.join('?' * len(assemblies))}) AND read_name = ? AND contig IS NOT NULL "
            f"AND start IS NOT NULL ORDER BY assembly, contig, start",
            list(assemblies) + [name]).fetchall()
    finally:
        con.close()
    out = []
    for asm in assemblies:
        mine = [r for r in rows if r[0] == asm]
        blocks = [r[1:-1] for r in mine if r[12] is not None and r[12] >= pid
                  and r[13] is not None and r[13] >= plen]
        best = {(b[1], b[2], b[3], b[10]) for b in _best_placements(blocks, tol)[0]}
        for r in mine:
            b = r[1:-1]
            ok = r[12] is not None and r[12] >= pid and r[13] is not None and r[13] >= plen
            out.append(dict(assembly=asm, contig=b[1], start=b[2], end=b[3],
                            strand="+" if b[4] == 1 else "-", pid=b[11], plen=r[13],
                            type=b[10], passes=bool(ok),
                            best=ok and (b[1], b[2], b[3], b[10]) in best))
    return out


def _best_by_read(blocks: list, tol: float) -> dict:
    """read -> (best placements, n) over one assembly's aln rows (`_best_placements`)."""
    by_read: dict = {}
    for b in blocks:
        by_read.setdefault(b[0], []).append(b)
    return {rn: _best_placements(bl, tol) for rn, bl in by_read.items()}


def _coverage_at_best(best: dict, segments: list, binsize: int, read_buckets) -> CoverageBins:
    """Marginal coverage: each read adds 1/n to the bins under each of its n best
    placements in this assembly (a supplementary piece of a best primary adds its
    primary's share), whichever of them lie in `segments`.

    Bins as `compute_coverage_bins` lays them out -- per segment, in (contig, start)
    order, the last of a segment cut short -- so bin i here is bin i of
    RegionData.bin_edges. Built from the aln rows already fetched, no BAM scan.
    """
    segments = sorted(segments, key=lambda s: (s.contig, s.start, s.end))
    bin_contig, bin_start, bin_end = [], [], []
    by_contig: dict = {}
    first = 0
    for seg in segments:
        n = math.ceil((seg.end - seg.start) / binsize)
        by_contig.setdefault(seg.contig, []).append((seg.start, seg.end, first))
        for j in range(n):
            gs = seg.start + j * binsize
            bin_contig.append(seg.contig)
            bin_start.append(gs)
            bin_end.append(min(seg.end, gs + binsize))
        first += n
    n_bins = len(bin_contig)
    name_to_bucket = {nm: i for i, bucket in enumerate(read_buckets) for nm in bucket}
    bins = [0.0] * n_bins
    bucket_bins = [[0.0] * n_bins for _ in read_buckets]
    bin_reads: list = [dict() for _ in range(n_bins)]
    read_names: list = []
    for rn, (kept, n) in best.items():
        if not n:
            continue
        w = 1.0 / n
        bi = name_to_bucket.get(rn)
        target = bucket_bins[bi] if bi is not None else bins
        ri = None
        for b in kept:
            for s0, s1, b0 in by_contig.get(b[1], ()):
                lo, hi = max(b[2], s0), min(b[3], s1)
                if lo >= hi:
                    continue
                if ri is None:
                    ri = len(read_names)
                    read_names.append(rn)
                for i in range(b0 + (lo - s0) // binsize, b0 + math.ceil((hi - s0) / binsize)):
                    target[i] += w
                    bin_reads[i][ri] = bin_reads[i].get(ri, 0.0) + w
    bin_offsets = np.zeros(n_bins + 1, dtype=np.int64)
    np.cumsum([len(s) for s in bin_reads], out=bin_offsets[1:])
    read_idx = np.fromiter((ri for s in bin_reads for ri in sorted(s)),
                           dtype=np.int32, count=int(bin_offsets[-1]))
    read_w = np.fromiter((s[ri] for s in bin_reads for ri in sorted(s)),
                         dtype=np.float64, count=int(bin_offsets[-1]))
    return CoverageBins(bins=bins, bucket_bins=bucket_bins, bin_contig=bin_contig,
                        bin_start=bin_start, bin_end=bin_end, bin_offsets=bin_offsets,
                        read_idx=read_idx, read_names=read_names, read_w=read_w)


def _placements(best: dict, segs_info: list, read_pos: dict) -> Placements:
    """`Placements` from `_best_by_read` output: the same reads and weights
    `_coverage_at_best` bins, kept as intervals instead of bins."""
    by_contig: dict = {}
    for k, (contig, s, e, off) in enumerate(segs_info):
        by_contig.setdefault(contig, []).append((s, e, off, k))
    rows = []
    for rn, (kept, n) in best.items():
        g = read_pos.get(rn)
        if not n or g is None:
            continue
        for b in kept:
            for s0, s1, off, k in by_contig.get(b[1], ()):
                lo, hi = max(b[2], s0), min(b[3], s1)
                if lo < hi:
                    rows.append((g, off + lo - s0, off + hi - s0, 1.0 / n, k))
    rows.sort()
    a = np.array(rows, dtype=float).reshape(-1, 5)
    read = a[:, 0].astype(np.int32)
    n_reads = len(read_pos)
    return Placements(p0=a[:, 1], p1=a[:, 2], w=a[:, 3], read=read, seg=a[:, 4].astype(np.int32),
                      first=np.searchsorted(read, np.arange(n_reads), "left").astype(np.int32),
                      last=np.searchsorted(read, np.arange(n_reads), "right").astype(np.int32))


MARGIN_BINS = 1000          # bins across the marginal's visible span


def nice_bin(span: float, target: int = MARGIN_BINS) -> int:
    """The bin size for a visible `span`: span/target rounded up to 1, 2 or 5 x 10^k
    bp, never below 1. Integer steps, so glyphs' REBIN_JS picks the same size."""
    raw = span / target
    if raw <= 1:
        return 1
    mag = 1
    while mag * 10 <= raw:
        mag *= 10
    return next(f * mag for f in (1, 2, 5, 10) if f * mag >= raw)


def coverage_gaps(segs_info: list, pl: Placements) -> tuple:
    """(left, right) on the axis of every stretch of a segment no placement covers:
    the coverage dropouts, exact to the bp whatever the bin size."""
    order = np.lexsort((pl.p0, pl.seg))
    left, right = [], []
    for k, (_, s, e, off) in enumerate(segs_info):
        reach = off                       # covered up to here
        for i in order[np.searchsorted(pl.seg[order], k, "left"):
                       np.searchsorted(pl.seg[order], k, "right")]:
            if pl.p0[i] > reach:
                left.append(reach)
                right.append(pl.p0[i])
            reach = max(reach, pl.p1[i])
        if off + (e - s) > reach:
            left.append(reach)
            right.append(off + (e - s))
    return np.array(left, float), np.array(right, float)


def marginal_bins(segs_info: list, pl: Placements, bands, n_bands: int,
                  a: float, b: float, binsize: int | None = None, flips=None) -> dict:
    """Bin `pl` over the axis range [a, b]: each row adds its weight to every bin it
    touches, in band `bands[read]`. Bins restart at each segment's start, the last of
    a segment cut short, as `_coverage_at_best` lays them out; only bins meeting
    [a, b] are made. The Python twin of glyphs' REBIN_JS, which redoes this in the
    page on every zoom; it draws the first frame.

    Returns left/right (axis), seg, contig, bin_start/bin_end (genome) per bin, and
    `series`: n_bands rows of per-bin values.
    """
    binsize = binsize or nice_bin(b - a)
    left, right, segk, contig, gs, ge, base = [], [], [], [], [], [], {}
    for k, (c, s, e, off) in enumerate(segs_info):
        n = math.ceil((e - s) / binsize)
        if off >= b or off + (e - s) <= a:
            continue
        j0, j1 = max(0, math.floor((a - off) / binsize)), min(n, math.ceil((b - off) / binsize))
        base[k] = (len(left) - j0, j0, j1)
        for j in range(j0, j1):
            l = j * binsize
            left.append(off + l)
            right.append(off + min(e - s, l + binsize))
            r = min(e - s, l + binsize)
            segk.append(k)
            contig.append(c)
            gs.append(e - r if flips and flips[k] else s + l)    # a flipped segment runs e -> s
            ge.append(e - l if flips and flips[k] else s + r)
    nseg = len(segs_info)
    b0 = np.zeros(nseg, np.int64)
    j0 = np.zeros(nseg, np.int64)
    j1 = np.zeros(nseg, np.int64)          # j0 == j1: no bins, nothing lands
    offs = np.array([si[3] for si in segs_info], float).reshape(-1)
    for k, (x0, x1, x2) in base.items():
        b0[k], j0[k], j1[k] = x0, x1, x2
    sg = np.asarray(pl.seg, np.int64)
    lo = np.maximum(j0[sg], np.floor((pl.p0 - offs[sg]) / binsize).astype(np.int64))
    hi = np.minimum(j1[sg], np.ceil((pl.p1 - offs[sg]) / binsize).astype(np.int64))
    keep = lo < hi
    band = np.asarray(bands, np.int64)[pl.read[keep]]
    diff = np.zeros((n_bands, len(left) + 1))
    np.add.at(diff, (band, (b0[sg] + lo)[keep]), pl.w[keep])
    np.add.at(diff, (band, (b0[sg] + hi)[keep]), -pl.w[keep])
    series = np.cumsum(diff, axis=1)[:, :-1]
    series[np.abs(series) < 1e-9] = 0.0      # what the running sums leave of an empty bin
    return dict(left=np.array(left, float), right=np.array(right, float),
                seg=np.array(segk, np.int32), contig=contig,
                bin_start=np.array(gs, np.int64), bin_end=np.array(ge, np.int64),
                series=series, binsize=binsize)


def _seg_at(offs: np.ndarray, p) -> np.ndarray:
    """Index of the segment holding axis position `p` (offsets ascending)."""
    return np.clip(np.searchsorted(offs, p, "right") - 1, 0, max(len(offs) - 1, 0))


def diagonal_order(rd: "RegionData", axis: str) -> tuple:
    """(order, flips) for the segments on `axis` ("x": asm2, "y": asm1), laid along
    the other axis so the reads run on one diagonal.

    A segment's place is the weighted median, along the other axis, of the reads'
    lines through it; it is flipped when those lines mostly run backwards. A line
    weighs its length over the read's line count, so a repeat read with many
    copies pulls no harder than a unique one. Partial rows do not vote. Segments
    no line reaches go last, in their current order.
    """
    L = rd.lines
    u0, u1, v0, v1 = (L.x0, L.x1, L.y0, L.y1) if axis == "x" else (L.y0, L.y1, L.x0, L.x1)
    segs = rd.segs_info(rd.asm2 if axis == "x" else rd.asm1)
    offs = np.array([x[3] for x in segs], float)
    use = np.asarray(L.layer) != LAYER_PARTIAL
    per_read = np.bincount(np.asarray(L.read_idx)[use], minlength=len(rd.read_names))
    k = _seg_at(offs, (u0 + u1) / 2)
    w = np.abs(u1 - u0) / np.maximum(per_read[L.read_idx], 1)
    v = (v0 + v1) / 2
    sign = np.sign((u1 - u0) * (v1 - v0))
    keys, flips = [], []
    for j in range(len(segs)):
        m = use & (k == j) & (w > 0)
        if not m.any():
            keys.append((1, j))
            flips.append(False)
            continue
        o = np.argsort(v[m])
        cw = np.cumsum(w[m][o])
        keys.append((0, float(v[m][o][np.searchsorted(cw, cw[-1] / 2)])))
        flips.append(bool(np.sum(w[m] * sign[m]) < 0))
    order = sorted(range(len(segs)), key=lambda j: keys[j])
    return order, flips


def arrange(rd: "RegionData", axis: str, order=None, flips=None) -> "RegionData":
    """`rd` with the segments on `axis` re-laid in `order` (default: `diagonal_order`),
    those flagged in `flips` (indexed by current segment) running backwards.

    Every axis coordinate moves with its segment -- lines, placements, indel marks
    -- so nothing is fetched again; genome coordinates are unchanged. Applies to a
    RegionData in its built layout (no flips yet on that axis).
    """
    if order is None:
        order, flips = diagonal_order(rd, axis)
    asm = rd.asm2 if axis == "x" else rd.asm1
    segs = rd.segs_info(asm)
    old_off = np.array([x[3] for x in segs], float)
    length = np.array([x[2] - x[1] for x in segs], float)
    flip = np.array(flips if flips is not None else [False] * len(segs), bool)
    new_off = np.zeros(len(segs))
    new_off[order] = np.concatenate([[0], np.cumsum(length[order])[:-1]]) if len(segs) else []
    pos_of = np.empty(len(segs), np.int32)
    pos_of[order] = np.arange(len(segs), dtype=np.int32)

    def move(p, k):
        rel = np.asarray(p, float) - old_off[k]
        return new_off[k] + np.where(flip[k], length[k] - rel, rel)

    L = rd.lines
    a0, a1 = (L.x0, L.x1) if axis == "x" else (L.y0, L.y1)
    k = _seg_at(old_off, (np.asarray(a0) + np.asarray(a1)) / 2)
    moved = dict(zip(("x0", "x1") if axis == "x" else ("y0", "y1"), (move(a0, k), move(a1, k))))
    lines = replace(L, **moved)

    pl = rd.placements[asm]
    sk = np.asarray(pl.seg, np.int64)
    q0, q1 = move(pl.p0, sk), move(pl.p1, sk)
    placements = {**rd.placements, asm: replace(pl, p0=np.minimum(q0, q1), p1=np.maximum(q0, q1),
                                                seg=pos_of[sk])}
    indels = dict(rd.indels)
    key = "x" if axis == "x" else "y"
    if len(indels.get(key, [])):
        pts = np.asarray(indels[key], float)
        indels[key] = move(pts, _seg_at(old_off, pts))

    new_segs = [(c, s0, e0, int(new_off[j])) for j in order for c, s0, e0, _ in [segs[j]]]
    field_name = "asm2_segs_info" if axis == "x" else "asm1_segs_info"
    return replace(rd, lines=lines, placements=placements, indels=indels,
                   flips={**rd.flips, asm: [bool(flip[j]) for j in order]},
                   **{field_name: new_segs})


def _compute_lines(best, read_names, read_pos,
                   asm1_segs_info, asm2_segs_info, asm1, asm2) -> ReadLines:
    """Thread each read between its best placements in asm1 and asm2 into plot-space
    line segments. `best[asm][read]` is `(placements, n)` from `_best_by_read`.

    Geometry helpers are reused verbatim from `qatool.plotters.reads_dotplot`; only
    the merging and matplotlib drawing are left behind.
    """
    y_lookup = _build_lookup(asm1_segs_info)   # asm1 -> Y
    x_lookup = _build_lookup(asm2_segs_info)   # asm2 -> X

    x0s, y0s, x1s, y1s, idxs, supps = [], [], [], [], [], []
    ct1s, ct2s, r1_0, r1_1, r2_0, r2_1, p1s, p2s = [], [], [], [], [], [], [], []
    n_dropped = 0
    for rn in read_names:
        a1 = best[asm1].get(rn, ([], 0))[0]
        a2 = best[asm2].get(rn, ([], 0))[0]
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
