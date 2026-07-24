 # Session Context — assemblyQA plotters.py

## Project overview

HiFi assembly QA tool. Given a candidate assembly region and a reference, it:
1. Recruits reads from a BAM-indexed SQLite DB (`SRRecruiter.py`)
2. Produces a `report_plot`: a combined k-mer dotplot + read-alignment overlay + marginal coverage plots

Entry point in notebooks: `testing_pipeline.ipynb`.

---

## Key files

| File | Role |
|------|------|
| `plotters.py` | All plotting logic (main work of these sessions) |
| `SRRecruiter.py` | Reads recruiter; calls `report_plot` |
| `SegmentStorage.py` | `Segment` dataclass + `SegmentCollection` (non-overlapping interval set) |
| `aln_utils.py` | CIGAR parsing (`cigar_stats`) |
| `testing_pipeline.ipynb` | Interactive test harness |

---

## plotters.py — function map

### `_build_axis_map(sc, assembly, binsize=None)`
Iterates `SegmentCollection`, returns `(segs_info, total_len)` where `segs_info` is
`List[(contig, start, end, offset)]`. When `binsize` is given, each segment's offset
contribution is `ceil(seg_len / binsize) * binsize` — so every segment boundary falls
exactly on a bin edge (no mixed bins).

### `alignment_plot(bamfile, segments, ..., binsize, orientation, return_max)`
Bar/area plot of read coverage along concatenated segments. Uses bin-aligned bin count
(`sum(ceil(seg_len/binsize) for seg in segments)`). Step is `'post'` everywhere (was
`'pre'` — caused a 1-bin left/down shift). Returns `(ax, max_coverage)` when
`return_max=True`.

### `dotplot_kmer_mosaic(seqs_a, seqs_b, k, binsize, ax, ...)`
**Mosaic** k-mer dotplot — NOT a single concatenated dotplot. Each pair `(seqs_a[i],
seqs_b[j])` is computed independently and placed in its cell of the combined matrix.
This eliminates mixed-bin artefacts at segment boundaries.

Algorithm:
- 2-bit encoding, forward + RC for seqs_a, forward-only for seqs_b
- Packed int64 `(hash << SHIFT | bin)` sorted arrays, `np.unique` to deduplicate
  `(hash, bin)` pairs (prevents tandem-repeat explosion)
- Repeat filter: k-mers with > `max_kmer_freq` occurrences in seqs_a are dropped
- Per-pair matching via `searchsorted`, ragged-range cumsum to expand a-side bins,
  `np.bincount` to fill sub-matrix
- `imshow` with `extent=(0, total_bins_b*binsize, 0, total_bins_a*binsize)` maps
  pixel coordinates back to base-pair space

### `reads_alignment_dotplot(sc, read_names, ..., binsize)`
Overlays read alignments as colored lines on top of the dotplot. Uses the same
bin-aligned `_build_axis_map` coordinate system so lines align with dotplot pixels.
Segment boundary lines: full-opacity for cross-contig, pale (alpha×0.5, lw×0.5) for
same-contig boundaries (recruiter always groups same-contig segments together).
Tick labels show only coordinates (`start:end`), no contig name.

### `report_plot(seq_asm1, seq_asm2, segments, read_names, db_path, ..., max_coverage=None)`
Main report function. Layout (bottom to top, left to right):

```
[ ax_top   ]   ← marginal coverage, asm2 (orientation="x"), size=2in
[ ax (dotplot) + ax_overlay (reads) ]
[ ax_xbar  ]   ← contig color bar (bottom), size=0.15in
[ ax_ybar  ]   ← contig color bar (left),   size=0.15in
[ ax_right ]   ← marginal coverage, asm1 (orientation="y"), size=2in
```

Steps:
1. `_build_axis_map` with binsize → bin-aligned offsets for both assemblies
2. Split concatenated sequences into per-segment lists (raw lengths, not padded)
3. `dotplot_kmer_mosaic` → k-mer background
4. Transparent `ax_overlay` inset → `reads_alignment_dotplot` on top
5. `make_axes_locatable` → colored contig bars (bottom/left); same-contig segments
   merged into one bar via `_contig_runs()`; caption inside bar if it fits, else legend
6. Marginal alignment plots (`ax_top`, `ax_right`); both coverage axes set to the
   same limit `m = max_coverage or max(maxx, maxy)` and the **same tick values**
   via `mticker.MaxNLocator(nbins=5, integer=True).tick_values(0, m)`

**`max_coverage` parameter**: pass the same value to multiple `report_plot` calls to
make marginal plots directly comparable across regions/assemblies.

---

## Coordinate system

All three layers (mosaic dotplot, reads overlay, marginal plots) share the same
coordinate space: `x ∈ [0, plot_w]`, `y ∈ [0, plot_h]` where
`plot_w = Σ ceil(seg_len_b / binsize) * binsize` and similarly for `plot_h`.
Every segment starts at an exact multiple of `binsize`. One pixel = one bin = `binsize` bp.

---

## Known past bugs (fixed)

| Bug | Fix |
|-----|-----|
| 30s slowdown with many matches | `np.unique` dedup of `(hash,bin)` pairs before accumulation |
| Mixed bins at segment boundaries | Replaced concatenated dotplot with mosaic (`dotplot_kmer_mosaic`) |
| Marginal plots 1-bin shifted | Changed `step='pre'` → `step='post'` in `alignment_plot` |
| Segment boundary lines off by 1px | Bin-aligned offsets: boundary = `segs[i+1].offset` exactly |
| Marginal coverage axes different ticks | Explicit shared `MaxNLocator` ticks applied to both axes |
| `coverage_analysis2` import error | Removed stale `from coverage_analysis2 import ref_segments` from `testing_pipeline.ipynb` |
| Last bin not drawn in marginal plots | `x = np.arange(len(bins)) * binsize` gave n x-values; `step='post'` needs n+1. Fixed: `x = np.arange(len(bins)+1) * binsize` + trailing `0.0` appended to every series so the last bin renders its full horizontal extent. |