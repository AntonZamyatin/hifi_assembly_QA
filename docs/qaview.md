# qaview — interactive design document

Status: Panel app implemented (§3–§11). **Being rebuilt** on a project file, a
FastAPI server and a CSS-grid page — see §13, written 2026-09-24, which supersedes
§3's choice of Panel and §7.5b's seed picker. Written 2026-07-27.

`qaview` is an interactive exploration frontend for the assembly-QA report figure.
It is a new package (`src/qaview/`) that reuses `qatool`'s data layer and adds a
browser-based renderer with hovering, highlighting, mouse zoom and live parameter
changes.

---

## 1. Scope and division of labour

The static `report_plot` is **not** being replaced. The two renderers split by purpose:

| | `qatool.plotters.report_plot` | `qaview` |
|---|---|---|
| Purpose | publication figures | interactive exploration |
| Output | PNG/PDF via matplotlib | live page via Bokeh/Panel |
| k-mer dotplot | always drawn (background) | off by default, on demand per viewport |
| Coordinate space | bin-padded (see §5) | exact bp |
| Driven by | function call with a `SegmentCollection` | seed segments + recruiter run |

Both sit on **one shared compute layer**, so numbers can never diverge between what
you explore and what you publish. The "export figure" button (§7.5) is literally a
call into `report_plot` with the current view state.

This resolves an earlier question of whether the k-mer dotplot is still needed now
that read-alignment dotplots exist. It is — but only in the publication path and in
zoomed-in on-demand views, not as an always-on background. See §8.

---

## 2. The sizing fact that drives the architecture

One region's view is small. A ~500 kb concatenated-segment span at ~25x with 20 kb
reads is order 10³ reads; coverage is ~10³ bins × 6 buckets. Serialized, a few
hundred KB.

**The whole view therefore fits in the browser.** Hover, highlight, crosshair and
zoom need zero server roundtrip. Only recruitment, refetching and reclassification
touch Python. Over a remote connection to the compute host this distinction is the
difference between fluid and unusable, and it rules out any stack whose interaction
model is "rerun the script".

---

## 3. Stack: Panel + Bokeh

Chosen because Bokeh maps onto this specific figure unusually directly:

- **`MultiLine` + `ColumnDataSource`** — one row per read, extra columns for name,
  pid, plen, class, contig coordinates. `HoverTool` surfaces any of them for free.
- **`selection_glyph` / `nonselection_glyph`** — selecting reads dims the rest with
  no custom code. This is the "highlight things of interest" requirement.
- **Shared `x_range` / `y_range`** — marginals and contig bars track the dotplot
  viewport automatically. This is what `make_axes_locatable` fakes statically today.
- **`WheelZoomTool`** — mouse zoom, built in.
- **`CrosshairTool` with a shared `Span` overlay** — cursor projection across linked
  figures (§7.1).
- **`image` glyph** — drops the on-demand k-mer mosaic straight in.

Panel adds the widget layer, and critically lets **the same code run both as a
notebook widget and as a served app** (`panel serve`, port-forwarded), so the
notebook-vs-webapp decision stays open.

Everything is Python; qaview imports `qatool` rather than reimplementing DB access
or classification.

### Alternatives considered

| Stack | Why not |
|---|---|
| Plotly / Dash | Per-read hover forces either thousands of traces or None-separated vertices that lose read identity; highlighting needs a restyle roundtrip. |
| FastAPI + deck.gl / canvas | Best ceiling, genome-browser-grade — but weeks of JS for a payload this small. Revisit only if many regions must render at once. |
| Streamlit | Rerun-everything model; every interaction re-hits multi-GB SQLite. |
| ipympl | Reuses all current code but has poor tooltips and slow remote re-render. |

### Dependencies

`src/qaview/` is auto-discovered by the existing
`[tool.setuptools.packages.find] where = ["src"]`. Add to `pyproject.toml`:

```toml
[project.optional-dependencies]
view = ["bokeh", "panel"]
```

so the core install stays lean. Current env: Python 3.12.5, numpy 2.4, pandas 3.0,
matplotlib. No viz frontend installed yet.

---

## 4. Package layout

```
src/qaview/
  region.py    RegionData — recruiter output -> plain arrays (incl. bin->read index)
  glyphs.py    Bokeh figure builders (dotplot, marginals, contig bars)
  state.py     ViewState + export-to-report_plot
  app.py       Panel layout, widgets, tiered callbacks
  __main__.py  `panel serve` / `python -m qaview` entry
```

### What is reused, extracted, rewritten

| Reuse as-is | Extract from rendering | Rewritten |
|---|---|---|
| `plotters/db.py` — all three fetchers | `alignment_plot` binning loop → `compute_coverage_bins()` (§6) | all matplotlib rendering |
| `plotters/classify.py` — both classifiers | `_build_axis_map` minus bin padding (§5) | `make_axes_locatable` layout |
| `plotters/reads_dotplot.py` geometry — `_read_span`, `_read_to_ref`, `_map_point` | | `report.py` orchestration |
| `plotters/colors.py`, `_coord_ticks` | | |
| `plotters/kmer.py::dotplot_kmer` (on-demand only) | | |
| `SRRecruiter.Recruiter.run` | | |

Only `alignment_plot` needs real surgery — it is the one place where BAM iteration
and `stackplot` are fused.

`_merge_segments` (`reads_dotplot.py:85`) is **not** reused. It merges colinear
blocks at `tol=binsize`, which is resolution-dependent and would need recomputing on
every zoom. It was only ever a static-rendering optimisation; a few thousand
unmerged segments is nothing for Bokeh. qaview draws unmerged segments and encodes
depth with alpha overplotting.

---

## 5. Coordinate system

qaview uses **exact bp offsets** — no bin padding.

The bin-aligned padding in `_build_axis_map` (`ceil(seg_len/binsize)*binsize`) exists
solely so k-mer mosaic pixels never straddle a segment boundary. With the mosaic out
of the default path, the padding buys nothing and costs a whole class of off-by-one
bugs (see the fixed-bug table in `plotters.md`).

Consequences:

- Segments are still binned **individually** for coverage (`ceil(seg_len/binsize)`
  bins each), but the final bin of each segment is clipped to the segment end.
  Bokeh draws coverage from explicit `(left, right, top)` quad arrays, so a
  variable-width final bin is trivial — unlike matplotlib's `stackplot`, which needs
  a uniform `x = arange(n+1) * binsize`.
- The on-demand k-mer layer is placed with an explicit `image` extent per viewport,
  so it needs no global bin alignment either.
- Read lines are bp-native already and unaffected.

### Coordinate inverse (new)

`axis.py` has `_bins_to_genome_intervals` but no scalar inverse. The crosshair
readout needs one:

```python
def plot_to_genome(segs_info, pos) -> tuple[str, int] | None:
    """Plot coordinate -> (contig, genomic position)."""
```

A crosshair reading `1,234,567` in concatenated space is useless; it must read
`SUPER_1:1,801,234`.

### ViewState stores genomic, not plot, coordinates

`ViewState` records `(contig, start, end)` tuples, thresholds, toggles and the zoom
window in **genome space**. This is what lets the two renderers each build their own
coordinate mapping — qaview exact-bp, `report_plot` bin-padded — from the same state.

---

## 6. Data layer: `RegionData`

Built once per recruitment result. Plain arrays and DataFrames, no matplotlib
artists, JSON-serialisable for the browser.

```
RegionData
  segments        SegmentCollection (recruiter output)
  asm1_segs_info  [(contig, start, end, offset)]  exact-bp
  asm2_segs_info  [(contig, start, end, offset)]
  reads           DataFrame: read_name, class, pid1, pid2, plen1, plen2, ...
  lines           arrays x0, y0, x1, y1, read_idx, is_supplementary
  coverage        per-assembly: bucket_bins [n_buckets][n_bins], bin_left, bin_right
  bin_reads       CSR index: bin_offsets[n_bins+1], read_idx[]   <-- see below
```

### `compute_coverage_bins()` must emit read identity

This is the single requirement that has to be built in from day one, because
retrofitting it changes the function's signature and memory profile.

`alignment_plot` currently accumulates counts only — `target[i] += w` at
`coverage.py:79-85` — and discards which reads landed in which bin. Feature §7.4
(hover a coverage bin, highlight its reads) needs that mapping.

Emit it CSR-style: `bin_offsets` (length `n_bins+1`) plus a flat `read_idx` array.
Measured on the hg002 region above (1980 bins across both assemblies at 35x): 116,862
int32 entries, ~456 KB to the browser — more than the ~100 KB first estimated, but
still comfortably client-side. Hover
then becomes a `CustomJS` that sets `source.selected.indices` on the read
`MultiLine` — instant, no roundtrip. The reverse direction (hover a read, light up
its bins) falls out of the same index for free.

### 6.1 Boundary clipping

Blocks overhanging a displayed segment are clipped to the visible part instead of
being discarded. The helpers live in `qatool.plotters.reads_dotplot` so both
renderers share one implementation:

| helper | role |
|---|---|
| `_ref_to_read` | inverse of `_read_to_ref` within a block |
| `_overlapping_segment` | segment with the greatest overlap (replaces the containment test) |
| `_visible_read_span` | read-coordinate interval whose reference image is on screen |
| `_map_in_segment` | reference position → plot coordinate, clamped |

**Clipping happens in read space, never per axis.** A line's two endpoints must refer
to the same stretch of read; clipping each axis independently would leave them
describing different spans and silently fake the slope. So each side contributes the
read interval it can show, and the intersection is what gets drawn:

```
ro0 = max(r0a, r0b, v1.lo, v2.lo)
ro1 = min(r1a, r1b, v1.hi, v2.hi)
```

`reads_overlay` needs the same intersection but has no read coordinates — its query
returns no clip fields. It therefore parameterises both alignments by `t ∈ [0, 1]`
along the read, which is the correspondence the overlay already assumes when it draws
a straight line between the two assemblies, and intersects the visible `t` ranges.

Result on the hg002 region: qaview 420 → **832 lines, 0 dropped, all 209 reads
represented**; the `report_plot` overlay 15 → **18** discriminating reads drawn.
`ReadLines.n_dropped` now counts only genuine misses (a block on no displayed segment,
or no shared visible read span), and the >25% warning still fires on those.

### Known data caveats to carry over

- **Reads dropped at segment boundaries — found, measured, fixed.** The original rule
  required a block sit entirely inside one displayed segment
  (`s <= ref_start and ref_end <= e`) and discarded it whole otherwise. For read
  length L over segment length S that loses `2L / (S + L)` of reads. On hg002
  `h2tg000020l:61,889,136-61,929,137` (4 recruited segments of ~41.5 kb) it cost
  **103 of 209 reads — 49.3%**, median overhang 7.1 kb, max 28 kb; below S = L it
  loses everything, so the 6.7 kb segment in the repo-root `report_plot.png` shows no
  reads at all. The loss was also biased, falling entirely on edge-spanning reads.
  This was **pre-existing, not a qaview regression** — `reads_overlay._map_pos`
  applied the identical rule and dropped exactly the same 103 reads.
  **Fixed in `qatool`** (§6.1) so both renderers recover them.
- **Coverage dropouts read as absence.** A read contributes a line only if it aligns
  to both assemblies; a dropout region renders blank, which looks like "no
  correspondence" rather than "no reads". The marginal coverage panels disambiguate,
  which is an argument for keeping them always visible.
- **Negative pid.** `cigar_stats` computes `pid = 1 - nm/q_aln_len` where NM counts
  deletions but `q_aln_len` excludes them, so deletion-heavy alignments go negative
  (min −1.80 observed in `mPleAur1.aln.sqlite`). Harmless while thresholding at
  ~0.99, but any ranking, averaging or colour-mapping on pid must clamp.
- **Secondary alignments dominate and are real.** In `mPleAur1.aln.sqlite`,
  secondaries are 108.7M of ~190M rows (57%) with mean pid 0.960; only 3% sit at
  pid=1.0, so NM tags are present rather than defaulting. Alignment uses
  `-N 20 -p 0.5 --secondary=yes` and `build_aln_db` stores every mapped record, so
  repeat and multi-locus structure is genuinely available to the read dotplot.

---

## 7. Feature designs

### 7.1 Cursor projection (ruler lines)

`CrosshairTool` sharing a `Span` overlay across the dotplot and both marginals
projects the cursor onto all linked panels. Paired with `plot_to_genome` (§5) for a
`contig:position` readout.

**Verified** against bokeh 3.9.2: one `CrosshairTool` instance carrying a
`(Span(dimension="width"), Span(dimension="height"))` overlay, added to the dotplot
and both marginals, gives a linked ruler. `tests/test_glyphs.py::test_shared_crosshair`
asserts the panels share one instance — equal-but-separate tools would silently
un-link the ruler.

### 7.2 Zooming

`WheelZoomTool` plus shared ranges. No server recomputation: read lines are vector
geometry, coverage is already binned, and merging is skipped (§4). The only
zoom-triggered server work is the optional k-mer layer (§8), which is explicit.

### 7.3 Interactive parameter change

Panel's `value_throttled` fires on slider release rather than during drag — matching
the decision to keep one Python implementation of classification rather than porting
it to JS. See §9 for the cost tiers.

### 7.4 Coverage-bin hover → read highlight

Backed by the CSR `bin_reads` index (§6). `CustomJS` on hover writes
`source.selected.indices`; `nonselection_glyph` dims everything else. Client-side and
instant.

Three implementation notes from building it:

- **An invisible per-bin hover target is required.** Hovering the stacked coverage
  quads directly reports whichever bucket lies under the cursor, and misses entirely
  where that bucket is zero. `coverage_figure` therefore appends one full-height
  transparent quad per bin (`hit_renderer`), which always resolves the bin and whose
  `selection_glyph` doubles as the bin highlight for the reverse direction.
- **Two CSR indexes, not one.** bin → reads comes from `compute_coverage_bins`;
  read → line rows comes from `ReadLines.by_read()`, since one read contributes
  several line segments (one per block pair).
- **The two indexes use different read numbering.** `CoverageBins.read_names` is
  per-assembly and holds only reads that reached a bin, while `RegionData.read_names`
  is the global sorted list the line rows index into. The callback builds a lookup
  table once and caches it on `window`. `test_coverage_read_names_are_a_subset`
  guards the assumption that every coverage read exists globally.

### 7.5 Publication export

**Do not export from Bokeh.** PNG export requires selenium + a chromedriver, and the
SVG backend mangles text — not worth fighting for a figure going into a paper.

Instead the export button serialises `ViewState` and calls `report_plot()`, which is
already debugged and produces vector/300 dpi output. Round-tripping through genomic
coordinates (§5) makes the two coordinate systems a non-issue.

`visible_segments()` clips the exact-bp axis map to the live Bokeh viewport and emits
genomic `Segment`s; `report_plot` then builds its own bin-padded map from those. Two
details that matter:

- **Round, do not truncate.** Bokeh range ends arrive as floats, so `int()` would
  shave a base off every segment on every export.
- **Export defaults `binsize=None`.** qaview auto-sizes to ~1000 bins with no k floor,
  so on a small region it lands *below* k (7 vs k=15 on the test fixture) and the
  k-mer mosaic would refuse to build. Letting `report_plot._resolve_binsize` choose
  applies its `max(k, …)` rule. The export panel still offers an override.

Rendering runs through `ExportJob` — the same `SubprocessJob` base as recruitment
(`qaview/jobs.py`) — so matplotlib never draws on a server thread and a slow k-mer
mosaic cannot block the page.

**Interpretation note:** exporting a zoomed sub-region changes multi-mapper weighting.
`compute_coverage_bins` splits a read's weight across its qualifying alignments *within
the displayed segments*, so restricting the view to one of two homologous segments
lifts those reads from 0.5 to 1.0 — on the hg002 region, marginal coverage goes from
~34x to ~60x. That is correct for a restricted view, not a bug, but the two figures are
not directly comparable.

---

### 7.5b Seeds from QUAST

The sidebar's first section lists QUAST misassemblies as seeds
(`qaview/seeds.py`), read from `<base>/quast/<run>/misassemblies.tsv` -- the table
the hg002_hic notebooks are generated from. Each row is named as the notebooks
name it, `<contig>_<type>_<i:03d>` with `i` counting rows of one (contig, type) in
file order, so a seed matches its notebook figure; the table shows a short form
(`h1tg000005l reloc1`) with a filter box.

- **breakpoint** -- `bp_start`; the full range is in the status line on picking.
- **ΔS nb** -- the notebook's pair-HMM score from
  `<base>/analysis/<run>/phmm_scores/<region>.json`, when present. It used the
  notebooks' recruitment parameters, which need not match the sidebar's. `reload`
  picks up scores written since the page opened.
- **✓** -- a recruitment for this seed with the current parameters is cached.

Picking a row fills contig/start/end (the +/-20 kb QUAST window). A cached seed
loads at once; otherwise press Recruit.

### 7.6 Pair-HMM score: log-ratio colouring and Delta S

**Layout.** The colour mode and the C slider sit above the plot, beside the status
line, so the colouring is never hidden in a collapsed section. The page is fitted
to the browser window: a small `Viewport` component reports the window size (and
debounced resizes), and `report_layout(max_width, max_height, margin_size)` fits
the dotplot into what is left, keeping square pixels. Marginals shrink from 200 to
120 px on windows under 1100 px tall. The sidebar and the main pane scroll on
their own; the page never does.

Each read's pair-HMM log-ratio between the two candidates (`qamodel.score`) is
drawn on the read itself, and the locus score sits above the plot.

**Scoring is its own tier.** About 30 s for the HG002 seed locus's 106 reads on
the compiled engine, so it runs like recruitment: a subprocess
(`qaview.score_worker`), cached as `score_<hash>.pkl`, started automatically when
a recruitment result arrives. The key is the recruitment key plus everything that
changes a ratio -- lambda's reference (asm2's primary reads over asm2's diploid
length), window sizing, band width, pair-HMM parameters.

**The clamp C is display-only.** `LocusScore.at(C)` re-evaluates Delta S from the
stored per-read ratios, so the C slider above the plot updates the panel and the
colour saturation instantly and never rescores.

**Colour: `log-ratio`.** A diverging scale on the palette's own hues -- asm1 in the
red of `asm1_only`, asm2 in the green of `asm2_only` -- through a neutral grey.
Symmetric-log rather than linear, `t = sign(r) log1p(min(|r|, C)) / log1p(C)`,
because at the seed most reads sit within +/-1 nat while the decisive ones reach
-27, and a linear 0..200 scale would render those almost neutral. Reads placeable
in only one candidate saturate. Unscored reads are a pale grey; dropped reads say
so in the hover. It is the default mode: until the score arrives the lines and
marginals draw the support classes, then switch once scoring finishes.

**Rows are never reordered.** The hover and selection callbacks address rows by
position, so decisive reads cannot simply be drawn last. Width (`1 + 2.5|t|`)
and opacity (`0.3 + 0.7|t|` of the base) carry the emphasis instead, which lets
indifferent reads recede. Hidden layers stay at alpha 0.

**Hover** gains a `log-ratio` line: the unclamped value in nats, `+/-inf` with
which candidate for one-sided reads, or why a read was not scored.

**Marginals follow the reads.** With colour set to `log-ratio` and a score in
hand, the coverage marginals stack by log-ratio band instead of support class:
indifferent (|r| < 1), weak (1-10 nats) and strong (>= 10 nats, about one Q40
substitution) for each candidate, and not scored. The bands are fixed in nats, not
scaled by C -- the clamp limits how much a read says, not which way it leans.
Before scoring finishes the marginals keep their class stacking.

The regrouping is exact and needs no BAM access. `CoverageBins` now stores each
read's weight in each bin (`read_w`, parallel to `read_idx`), and
`CoverageBins.rebucket` sums it under any grouping. Membership alone would not
do: a read with n tied alignments weighs 1/n per placement, and 2/n where two of
them fall in one bin -- which happens precisely in repeats. The log-ratio bands
sum to the same coverage as the class bands on both assemblies.

The stacking -- keys, colours, per-bin series, and the band of every read -- is
one `Stacking` object with two constructors, so the hover narrowing works for
either: hovering the "curated >= 10 nats" band isolates exactly the reads counted
in it.

Two readings to keep straight on the HG002 seed locus:

- *A strong band marks decisive reads, not the difference itself.* Bands are per
  read, and a read's coverage runs its whole ~15 kb, so the 7 bp microsatellite
  difference shows as a ~28 kb band around it. It also appears on the other
  haplotype, because the coverage track splits a read across every placement
  within `best_pid_delta` of its best, and 7 bp in 15 kb moves identity by only
  ~0.0005. At identity resolution the placements tie; the pair HMM separates them
  by 27 nats.
- *The pale "not scored" band is the containment rule.* Coverage counts every
  read overlapping a segment; scoring uses only reads inside one. Unscored
  coverage falls from 81% within 5 kb of a segment end to 7.5% in the middle --
  Proposition 8's "a third of the data at a 45 kb locus", made visible.

## 8. The k-mer layer

Off by default. Read lines are the interactive background because they are vector
geometry — free to zoom, nothing to recompute.

What the k-mer layer uniquely provides, and why it is kept on demand:

1. **The only read-independent evidence in the figure.** Every other layer derives
   from the same alignment DB, so a systematic mismapping in a repeat produces a
   confident wrong diagonal with nothing to contradict it. The k-mer background is
   what says "these two loci genuinely are the same sequence".
2. **Signal where coverage is zero** — distinguishes "no reads" from "no
   correspondence".
3. **Sub-block structure** — `_read_to_ref` interpolates linearly within a block, so
   indels inside a CIGAR are invisible in the read dotplot; the k-mer plot shows them
   as diagonal steps.

All three matter precisely when zoomed in, which is when computing them is cheap.
On demand, for one viewport, it is a plain `dotplot_kmer` call on two subsequences —
the mosaic machinery, `SHIFT` bit-packing and global bin alignment are all unneeded
because there is one segment pair in view.

The always-on mosaic stays in the `report_plot` publication path unchanged.

---

## 9. Parameter tiers

| Tier | Parameters | Cost | UI treatment |
|---|---|---|---|
| **Recruitment** | `merge_gap`, `extension_window`, `extension_k`, `ext_max`, `ext_rounds`, recruiter `pid`/`plen` | minutes | explicit "Recruit" button + progress |
| **Fetch / classify** | display `pid`/`plen`, `best_pid_delta`, `pid_tie_eps`, support vs occurrence | ~0.2–2 s | `value_throttled` |
| **Display** | layer toggles, colours, zoom, selection, binsize for display | instant | client-side |

Getting this boundary right is what determines whether the tool feels responsive.

**Measured on hg002 `h2tg000020l` (4 segments, 209 reads):** Tier 2 rebuilds a whole
`RegionData` — BAM read collection, two DB fetches, classification and both coverage
bins — in **0.7–1.0 s**, and Tier 3 redraws in **0.33 s**. A permissive-superset cache
was therefore not needed; `value_throttled` straight onto `RegionData.build` is fast
enough, and keeping one code path avoids a second, divergent filter implementation.
Revisit if a region an order of magnitude larger proves sluggish.

Tier 3 must never rebuild `RegionData` — `tests/test_app.py::test_tier3_does_not_rebuild_region`
asserts it, since a stray rebuild would still be *correct*, just 3x slower, and so
would never be noticed by eye.

---

## 10. Recruiter integration

`Recruiter.run(seed) -> (SegmentCollection, Set[str], DataFrame)` is the entry point.
The user supplies seed segments; qaview runs recruitment and visualises the full
recruited picture.

- **Subprocess, not a thread.** `qaview/recruit.py` runs recruitment via
  `python -m qaview.recruit_worker`. A thread would have needed
  `check_same_thread=False` for the SQLite connection opened at `SRRecruiter.py:45`,
  plus a global `sys.stdout` redirect to capture progress, and would not have been
  killable. A subprocess gives all three for free: the connection is opened and used
  in one process, progress is just the child's stdout, and Cancel is `terminate()`.
- **Caching.** The worker pickles `{"segments", "reads"}` to a path keyed by a SHA1 of
  the seed plus every recruitment parameter (`RECRUIT_KEYS`). Display and
  classification parameters are deliberately excluded from the key, so Tier-2/3
  changes can never invalidate it — asserted by
  `test_recruit_key_ignores_display_params`.
- **Progress.** The child inherits `PYTHONUNBUFFERED=1`; the parent pumps its stdout
  into a bounded deque that a 500 ms Panel periodic callback drains into the log pane.
  Polling rather than pushing, because a subprocess has no in-process callback to hook.
- Round scrubbing (`_write_round_plot`, `max_rounds`, `ext_rounds`) remains available
  for later; not implemented.

---

## 11. Phasing

| Phase | Deliverable |
|---|---|
| **0** | **DONE** — Compute layer: `RegionData`, `compute_coverage_bins()` with read identity, `plot_to_genome`. Verified structurally against synthetic fixtures (`tests/test_compute_layer.py`); a regression oracle holding the pre-refactor binning asserts `alignment_plot` is unchanged. Not yet run against real data. |
| **1** | **DONE** — Static Bokeh report in `qaview/glyphs.py`: read lines + marginals + contig bars, linked ranges, wheel zoom, shared crosshair. Rendered via `scripts/qaview_render.py`; structure covered by `tests/test_glyphs.py`. Still synthetic-data only. |
| **2** | **DONE** — Hover and selection: read tooltips, tap-to-highlight, coverage-bin → read highlight and the reverse, live coordinate readout (§7.1, §7.4). All client-side. `tests/test_interaction.py`. |
| **3** | **DONE** — `qaview/app.py` (Panel, three tiers), `qaview/recruit.py` + `recruit_worker.py` (cached, cancellable subprocess job), `qaview/__main__.py`. `tests/test_app.py`. |
| **4** | **DONE** — `qaview/state.py` (`ViewState`, `visible_segments`, `ExportJob`) + `export_worker.py`, wired into the sidebar. `tests/test_export.py`. |
| later | On-demand k-mer layer, recruitment-round scrubbing. |

**Phase 1's layout risk is retired.** `gridplot` reproduces the `make_axes_locatable`
proportions with explicit pixel sizes (`MAIN_WIDTH=1100`, `MARGIN_SIZE=200`,
`BAR_SIZE=15` at 100 px/in), and `main_height()` scales the dotplot by
`plot_h / plot_w` so pixels stay square — without it, diagonals sit off 45° on any
region whose two assemblies differ in length.  Two findings worth recording:

- Panels must share the **same** `Range1d` objects, not equal copies, or zoom
  un-links.  Asserted by identity in `test_ranges_are_linked`.
- `merge_tools=True` makes bokeh pick one active drag/scroll tool across panels that
  have different ones — it warns, and silently changes which panel the wheel drives.
  The layout therefore keeps a toolbar on the dotplot only.

---

## 12. Open questions

- Notebook widget vs served app as the primary mode — deferred by choosing Panel;
  revisit once Phase 3 exists.
- Whether display `binsize` should be a user control or derived from zoom level.
- ~~Whether the segment-boundary drop rate warrants clipping partially-contained
  blocks.~~ **Resolved: clipped in `qatool`** — see §6.1.
- **`pid_tie_eps` must be strictly less than `best_pid_delta`.** In
  `classify.py:167-175` the tie branch requires `abs(d) <= best_pid_delta` while a
  best call inside it requires `abs(d) > pid_tie_eps`; equal values make that a
  contradiction, so every non-common read falls through to `other` and the classifier
  can only ever emit common/other. On the hg002 region at the study's
  `pid_tie_eps = best_pid_delta = 0.001` this produced **zero** discriminating reads;
  at 0.0001 it gives 18, at 0.00001 it gives 30 (all the signal lives in
  `0 < |d| <= 0.001`). Worth a guard/warning in the classifier — not yet added.

---

## 13. Rebuild: project file, FastAPI server, CSS-grid page

Written 2026-09-24. Replaces Panel as the page framework; Bokeh stays as the plot
engine and `glyphs.py`, `region.py`, `score.py`, `seeds.py`, `jobs.py` are reused.

### 13.1 Why

Every layout problem in the Panel app came from Panel, not Bokeh: the window had
to be measured through a custom JS component, chrome and margins were guessed pixel
constants (`CHROME_HEIGHT = 280 + BAR_SIZE`), forms overflowed the sidebar, controls
hid in collapsed accordion sections, and Panel's shadow DOM defeats headless
screenshots and clicks. All plot interaction (hover, cross-highlighting, crosshair,
readout) is already client-side `CustomJS`, so it needs no Bokeh server and moves
unchanged. The workflow also changed: the batch showed that one QUAST run is ~250
seeds collapsing to ~100 loci, so the app should work through a whole comparison,
not one hand-typed seed at a time.

### 13.2 The project file

One TOML file describes one comparison of two candidates. Paths are relative to the
file. Opening it is the only thing the user does; everything else is derived.

```toml
name = "hg002_hic"

[asm1]                                    # Delta S > 0 favours asm1
name  = "hifiasm_hic"
fasta = "../bgz/hg002.hifiasm_hic.fasta.bgz"
bam   = "../alignments/hg002.hifiasm_hic.sorted.bam"

[asm2]
name  = "curated"
fasta = "../bgz/hg002.curated.fasta.bgz"
bam   = "../alignments/hg002.curated.sorted.bam"

[data]
db = "../dbs/hg002.aln.sqlite"
library_assembly = "curated"              # lambda reference; fixed across loci

[[seeds]]                                 # any number of sources
format   = "quast"                        # misassemblies.tsv from extract_misassemblies.py
path     = "../quast/hg002_hic/misassemblies.tsv"
assembly = "hifiasm_hic"                  # the candidate the coordinates are on

[[seeds]]
format   = "bed"                          # contig, start, end[, name]; 0-based half-open
path     = "extra_seeds.bed"
assembly = "hifiasm_hic"

[recruitment]                             # optional; defaults are the notebooks' values
pid_threshold = 0.997
plen_threshold = 0.999
best_pid_delta = 0.001
merge_gap = 50000
extension_window = 100
ext_max = 300000
ext_rounds = 100
max_rounds = 100

[resources]
recruit_workers = 32                      # seeds recruited at once
score_workers = 32                        # processes scoring one locus

[output]
dir = "../analysis/hg002_hic"             # recruit/, loci/ (score cache), state
```

**Seeds are a source-neutral list.** Each `[[seeds]]` entry is read by a loader for
its `format` into `Seed(name, source, assembly, contig, start, end, meta)`. QUAST
names keep the notebooks' `<contig>_<type>_<i:03d>` so rows match figures and
batch scores; BED names come from column 4 or `<contig>_<start>_<end>`. `meta`
carries source-specific fields (type, breakpoint, inconsistency) for display only.
New sources are a new loader, nothing else.

**Recruitment parameters belong to the project.** Every seed is recruited in the
background with them, so they are shown read-only in the page; changing them means
editing the file and reopening. `prepare_data.sh` will later write this file,
since it knows everything except the seeds.

### 13.3 Opening a project

1. Validate: every file exists, both assembly names occur in the database.
2. Load all seed sources; seeds with identical (assembly, contig, start, end) from
   different sources are one seed with several names.
3. **Recruitment queue.** Each seed not already cached is recruited in a worker
   process (`recruit_workers` at a time; seeds are independent). Results go to
   `output/recruit/recruit_<key>.pkl`, keyed by the seed, the data files and the
   recruitment parameters, so reopening reuses them.
4. **Loci.** A finished seed gets `locus_id(segments, reads)`; seeds with the same
   id are one locus, one table row, listing its seeds. Rows appear as seeds finish.
5. **Scoring queue.** Each distinct locus is scored once, automatically, through
   `phmm_likelihood_report(cache_dir=output/loci)` -- the cache the notebooks use,
   so loci the batch already scored arrive with Delta S at once. Scoring runs in a
   subprocess (it forks `score_workers` of its own; forking from the threaded
   server process is avoided). A locus the user opens moves to the front.
6. State is only ever derived from these files, so a crash or restart loses at
   most the jobs in flight.

### 13.4 Page

```
+----------------+-------------------------------------------+---------------+
| project        |  locus . Delta S . status                 |  (right       |
|  name, queue   +-------------------------------------------+   sidebar,    |
|  progress      |                                           |   reserved)   |
|                |                                           |               |
| loci table     |    dotplot + marginals + contig bars      |               |
|  Delta S,      |    (fills the remaining space)            |               |
|  reads, span,  |                                           |               |
|  seeds, status |                                           |               |
|                |                                           |               |
| plotting       |                                           |               |
|  display pid,  |                                           |               |
|  plen, deltas, |                                           |               |
|  binsize,      |                                           |               |
|  colour-by,    |                                           |               |
|  clamp C,      |                                           |               |
|  layers        |                                           |               |
|                |                                           |               |
| recruitment    |                                           |               |
|  (read-only)   |                                           |               |
+----------------+-------------------------------------------+---------------+
```

- CSS grid, three columns. Both sidebars scroll on their own and can be collapsed;
  the plot column takes what is left. One locus at a time.
- The **right sidebar** is reserved for per-locus and per-read information; its
  contents are to be designed. It ships as an empty, collapsible panel.
- **Sizing:** the page measures the plot cell (ResizeObserver) and asks the server
  for figures at that size; `fit_main` keeps square pixels server-side as now.
  Resizing re-requests after a 300 ms pause.
- **Plotting controls** are the Panel app's tiers 2 and 3: display pid, plen,
  best_pid_delta, pid_tie_eps (with the tie warning), binsize, colour-by, clamp C,
  layers. Tier 2 rebuilds `RegionData` server-side (~1 s, cached per locus and
  settings); tier 3 re-renders only.

### 13.5 Server

FastAPI + uvicorn, one process, reached by ssh port-forwarding.

| endpoint | returns |
|---|---|
| `GET /` | the page (static HTML, JS module, CSS) |
| `GET /bokeh/...` | BokehJS from the installed `bokeh` package -- no internet needed |
| `GET /api/project` | names, files, seed sources, recruitment parameters, resources |
| `GET /api/loci` | the table: id, seeds, segments, span per candidate, reads, status, Delta S |
| `GET /api/events` | server-sent events: queue progress and row updates |
| `GET /api/loci/{id}/view?w=&h=&...` | Bokeh `json_item` of the report layout for that size and those settings, plus the score summary |

Launch: `venv/bin/python -m qaview <project.toml> [--port 5006]`.

### 13.6 Reuse and removal

| module | fate |
|---|---|
| `glyphs.py` | reused; `report_layout` already takes a size box |
| `region.py`, `state.py` | reused |
| `recruit.py` | `run_recruitment` reused by the queue; cache key reused |
| `score.py`, `score_worker.py` | replaced by `phmm_likelihood_report(cache_dir=...)` in a subprocess |
| `seeds.py` | becomes the QUAST loader behind the seed interface |
| `app.py` | replaced by `project.py`, `pipeline.py`, `server.py` and `static/` |

### 13.7 Open

- Right-sidebar contents.
- Whether the read inspector (clicked read: per-read log-ratio, windows in both
  candidates, placements) goes there.
- Loci larger than ~5 Mb overplot at full view: aggregate reads at full zoom and draw
  individual reads only when zoomed in.
