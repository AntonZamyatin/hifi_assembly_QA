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

### 10.1 Extension strategies (`[recruitment] extension`)

Each round, before reads are collected, every segment's two edges are moved out
along its contig. Two strategies, chosen per project or per seed:

- **`legacy`** (default): crawl 100 bp windows (`extension_window`). A window with
  more than `extension_k` = 5 common reads, or with no uncommon read, seals the
  edge after snapping out to the nearest common read crossing it; otherwise the
  edge moves to the outermost start (end) of any qualifying read in the window.
  Seals are kept per exact segment, so a segment that changes is unsealed.
- **`anchor`** (experimental, branch `recruit-anchor`): edges move to anchor
  regions, so that a locus is **closed** -- every read touching it either lies
  inside its segments or overlaps them only on sequence both candidates share.
  - An **anchor region** is a best alignment of a *common* read (best pids tie
    within `best_pid_delta`, and the reference under its best placements is the
    same multiset of sequences in both candidates, either orientation) that no
    qualifying alignment **bridges** -- starts before it and ends after it. Every
    best locus of the read, in both candidates, must be unbridged; the read's
    best loci are the anchor's **copies**. A repeat read with three best loci in
    each candidate is an anchor with six copies.
  - An **anchored segment** ends at a contig end or at an anchor's outer end, the
    anchor inside it. A read crossing that edge from outside cannot bridge the
    anchor, so it ends within it: it overlaps the segment only on the anchor's
    sequence, the same in both candidates. Such reads are left out (no boundary
    reads any more); a read starting in the anchor and running inward lies inside
    the segment and is an ordinary locus read.
  - Each open edge moves to the **nearest** anchor region outward -- for a left
    edge the one with the largest start below it, for a right edge the smallest
    end above it -- or to the contig end. Coverage gaps hold no anchor and are
    crossed. `ext_max` limits one move; an edge stopped by it moves on next round.
  - Every other copy of that anchor becomes a segment of its own, sealed on the
    same outer side (the other side for a copy in reverse complement), and its
    inner edge moves the same way. This builds the other candidate's matching
    segment directly, and every copy of a repeat the anchor lies in.
  - Seals are positions, (assembly, contig, side), so they survive merges and
    rounds. The phase repeats until no edge is open, each edge position moving
    once per phase.
  - Recruited reads add every qualifying alignment as a segment, common reads
    included, as in legacy: HiFi-quality alignments are what tie related
    segments together, and the scorer must see each read's placements in both
    candidates. (Version 3 added none for common reads: on
    ptg000022l_28073999_28773000, 264 common reads whose curated copies were
    left outside the locus scored +C, and Delta S read +53.3k nats. Version 4
    added only their best placements: +721. Version 5, every alignment: 2.54 /
    2.52 Mb, 3,750 reads, 7 at +C, Delta S +785.)
  - `Recruiter.edge_log` records each edge decision (round, from, to, why: anchor,
    copy, contig_end, limit; and the anchor read) to compare strategies; it is
    saved as `edge_log` in the recruitment pickle. `extension_window` and
    `extension_k` are not used.
  - Anchor recruitments carry `extension_version` in their cache key
    (`pipeline.ANCHOR_VERSION`, 5 since common reads bring every qualifying alignment), so a change to this strategy
    redoes them without touching legacy ones.

The legacy strategy's cache key is the one it had before the option existed, so
existing recruitments stay valid; an anchor recruitment has a key of its own.

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
extension = "legacy"                      # or "anchor" (§10.1)

[resources]
recruit_workers = 32                      # seeds recruited at once
score_workers = 32                        # processes scoring one locus
synteny_threads = 32                      # minimap2 threads per genome alignment

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
4. **Loci.** Seeds that recruit the same read set are one locus, one table row,
   listing its seeds, even when their segments' ends differ (a seed's own
   coordinates, or an edge stopped by `ext_max`, can move an end with no read past
   it). The row keeps the segments of the seed that made it -- the first in project
   order when reopened -- and is named `locus_id(segments, reads)` from those. Rows
   appear as seeds finish. On HG002 (459 seeds) this makes 333 rows, where matching
   segments too made 352.
5. **Scoring queue.** Each distinct locus is scored once, automatically, through
   `phmm_likelihood_report(cache_dir=output/loci)` -- the cache the notebooks use,
   so loci the batch already scored arrive with Delta S at once. Scoring runs in a
   subprocess (it forks `score_workers` of its own; forking from the threaded
   server process is avoided). A locus the user opens moves to the front.
6. State is only ever derived from these files, so a crash or restart loses at
   most the jobs in flight.
7. **Re-recruitment.** A seed can be recruited again with parameters of its own
   (the Recruitment panel, `POST /api/recruit`). Only the values that differ from
   the project file are kept, per seed name, in
   `output/qaview_recruitment.json`; the seed's cache key uses its own values, so
   a reopened project recruits it the same way, and going back to the project
   file's values finds the old recruitment on disk at once. Requests run in a
   pool of their own (up to 4 workers), not behind the startup queue. The seed
   keeps its locus while it recruits; the result then moves it to the locus it
   makes (scored next), and a locus left with no seed is removed from the table.
   A failed request leaves the seed where it was, with the error on its row.

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
|  parameters,   |                                           |               |
|  Re-recruit    |                                           |               |
+----------------+-------------------------------------------+---------------+
```

- CSS grid, three columns. Both sidebars scroll on their own and can be collapsed;
  the plot column takes what is left. One locus at a time.
- The **right sidebar** holds, in its middle, the **alignment** of the stretches
  under a rectangle shift+dragged on the dotplot (a shift+click takes the single
  bins under the pointer); above and below are reserved. Each side is snapped
  outward to whole bins and kept to one segment, at most 200 kb a side. Both
  whole stretches are shown in one horizontally scrolling strip, hifiasm on top,
  laid out along their true correspondence: the asm1 stretch is located on asm2
  with edlib (HW, both orientations, a window centred on the rectangle), the
  path is polished with affine gaps (parasail) so an indel is one gap, and
  sequence of one stretch whose counterpart lies outside the other is shown
  hatched and not counted. `src/qaview/align.py`.
- **Show reads** (a toggle at the bottom of the right sidebar) adds the reads
  on each stretch -- above the top row, below the bottom row -- each read
  *projected* onto the alignment's columns through its CIGAR, as a pairwise
  alignment extends to a multiple one: aligned bases in their row base's column
  (bar, or dark where mismatched); row bases the read skips as a dark line;
  and where the row has gap columns (bases only the other assembly has), a read
  inserting there fills them as whole nucleotides (purple, compared with the
  other row), while a read inserting nothing shows a thin line, lacking them as
  its row does. Inserted bases with no gap columns to go in are an orange tick
  -- not necessarily bases neither assembly has: in a repeat the assemblies
  size differently, a read of a third length shows a tick against one row and
  a deletion against the other. In a tandem repeat an indel can sit anywhere
  along the repeat at the same score, so each read indel is walked through its
  equivalent positions in the row as shown (never past a mismatch or another
  indel, which would change the alignment) and placed on the alignment's gap
  columns when it can be, else at the leftmost, where the alignment's own gaps
  are left-aligned. The placements are exactly the ones the marginals count
  (every read at its best placements, display pid/plen, best_pid_delta),
  coloured and banded as the marginals are for the current colour mode, with
  band 0 next to the sequence on both sides. BAM records here carry no MD tag
  and secondaries no SEQ, so a read's sequence comes from its primary record;
  projected on its own reference a record's differences sum to NM exactly. The
  page measures the text's column grid (a run of letters placed in it) rather
  than deriving it from padding, so bars sit on the letters to a fraction of a
  pixel. `src/qaview/pileup.py`.
- **Remember** (beside Show reads, with an optional note) saves what the panel
  shows -- both stretches in genome coordinates, the settings, the alignment,
  and the reads if shown -- to `<output>/remembered/<seed>_<time>.json`, with a
  `.txt` of the alignment in 60-column blocks to read without the app. The
  snapshot is marked on the dotplot as a shaded yellow rectangle on any locus
  view whose segments contain both stretches (so it survives re-recruitment and
  shows for every seed recruiting the locus); hovering shows its note, and a
  click (not the end of a pan) restores it exactly as saved. With a remembered
  alignment open the button reads **Forget**, which deletes its files.
  The rectangles are BokehJS BoxAnnotations added by the page, changed in place
  rather than rebuilt, so remembering never re-renders the plot.
  `src/qaview/remember.py`.
- **Sizing:** the page measures the plot cell (ResizeObserver) and asks the server
  for figures at that size; `fit_main` keeps square pixels server-side as now.
  Resizing re-requests after a 300 ms pause. The grid's one row is the window's
  height (`minmax(0, 1fr)`), never its content's, so the sidebars scroll on
  their own and the plot cell is never taller than the window. Sidebars take at
  most 40% (left) and 30% (right) of the window's width, whatever was saved, so
  a narrow window shrinks them before the plot. The figure shrinks to any cell
  down to 160 px (smaller marginals below 450 px, and fewer ticks on them); the
  coordinate and status lines above it stay one line each. Only below Bokeh's
  own minimum sizes is the figure larger than its cell, and then the cell
  scrolls. With reads shown, the alignment strip is capped at 60% of the
  window's height and scrolls both ways, opening on the sequence rows.
- **Adding a seed.** "+" beside the table's filter opens a form: an assembly, a
  region as `contig:start-end` (1-based and inclusive, as the page shows
  coordinates; commas allowed) and an optional name (default
  `<contig>_<start>_<end>`, 0-based like BED). It is checked against the FASTA
  (contig, length) and the project's seeds (same region, name taken), saved to
  `<output>/qaview_seeds.json` -- not the project file -- which `load_project`
  reads as one more seed source every time, and recruited at once with the
  project's parameters; its locus is scored next and opens when it appears.
  `POST /api/seeds`.
- **Recruitment** (left sidebar) shows the open locus's seeds, each with its
  interval, and the parameters they were recruited with, editable; a value unlike
  the project file's is highlighted, with the file's value in its tooltip.
  **Re-recruit** recruits the ticked seeds (all, by default) again with the
  values in the form; **Project values** fills in the file's. While a seed
  recruits its row reads "recruiting"; when it is done the panel says where it
  went -- the same locus, or another one, with a link to it, and if the open
  locus lost its last seed the page opens the new one. A seed with parameters of
  its own is marked `*` in the table. The plotting controls' display pid, plen
  and best_pid_delta are separate from these: they filter what is drawn, not
  what is recruited.
- **k-mer similarity background** (`src/qaview/similarity.py`). For every pair
  of bins -- asm1 bin a (y), asm2 bin b (x) -- the share of a's distinct canonical
  21-mers also found in b estimates their identity, p_hat = share^(1/21) (a
  substitution kills about 21 k-mers). The share is taken both ways and the
  smaller counts (min-containment): one way only, a bin of few distinct k-mers --
  a uniform stretch of a centromeric satellite -- is contained in every bin
  holding its unit and lit a whole row (ptg000114l_translocation_001: 209 of
  1,000 rows lit across more than 40% of the columns; none now, the diagonals
  kept). Jaccard was considered: on equal-size bins it gives the same identity,
  but its union would count a cell's two-bin window, halving a split diagonal.
  Cells are drawn from the locus's recruitment pid up -- its own if its seeds
  were re-recruited with one, else the project's -- light at that floor, dark at
  100%: only sequence close enough for the recruitment to follow reads across,
  the diagonal and every repeat copy that near. The display pid does not move it:
  lowering that to see more reads used to lower the background's floor too, and
  flooded it. Bins are a regular grid along each
  axis (about 1,000 on the longer one), so the background is one image and costs
  the same to draw however many cells light up. The grids of the two candidates
  do not line up, so a cell is scored against the other axis's bin plus its
  better neighbour, in both directions. No k-mer is dropped for being frequent: past 200 million bin-pair
  counts, k-mers are hash-sampled (FracMinHash), which keeps the share unbiased.
  Computed in the background when a locus is first viewed (0.3 s on a 0.27 Mb
  locus, 9 s and 1.2 GB on the largest, 29 Mb), cached in `<output>/similarity/`;
  the page redraws when it arrives. In a satellite array every bin holds the same
  k-mer set, so the whole block is uniformly near 100%: true at this bin size,
  but without structure. (A minimap2 alignment of the segments was tried first:
  whole segments keep almost no repeat copies; 50 kb query windows keep them,
  13 s on 32 threads for 4.5 Mb, but the picture was less clear.)
- **Reads on hover.** "reads: on hover" (the default) draws no read line until
  a marginal bin is hovered, then only that bin's reads, at all their
  placements, a little thicker; "all" draws every read line as before. A frame
  then costs the background alone: on a 15.6 Mb locus 104-163 ms against
  240-485 ms with all 180,000 read lines (headless Chromium, canvas). The read
  lines are still sent with the view; fetching a bin's reads on hover instead is
  the next step for opening large loci faster.
- **Reads at their best placements.** Every read, whatever its class, is shown
  at its best placements only: in each assembly, its non-supplementary
  alignments passing the display pid/plen and within best_pid_delta of its best
  pid there -- over the whole assembly, not only the locus -- plus the
  supplementary pieces of a best primary (`region._best_placements`). The
  dotplot runs a read between its best placements in the two assemblies (a
  repeat read with k1 and k2 best copies draws k1 x k2 lines, its ambiguity);
  each marginal splits the read's weight of 1 over its n best placements in that
  assembly, a supplementary piece adding its primary's share, so a best copy
  outside the locus still takes 1/n and coverage reads as expected depth: a
  collapse doubles, a false duplication halves. The pileup under a rectangle
  alignment follows the same rule. Coverage is built from the aln rows already
  fetched, not a BAM scan. (Before, only reads voting for one assembly were drawn
  at their best placements, within pid_tie_eps; the others drew every pair of
  qualifying alignments, and n counted only the copies inside the locus.) On a
  6,000-read HG002 locus: 70,343 lines became 14,575, and building the region
  3.1 s became 1.7 s. pid_tie_eps now only splits asm1_best from asm2_best.
- **Partial alignments of one-sided reads.** Split and clipped alignments take
  no part in recruitment (supplementary records never recruit a read nor become
  segments: letting their pieces do so tied in unrelated regions through random
  partial alignments), nor in the default view. A read that only one candidate
  holds whole (asm1_only / asm2_only) may be a breakpoint read, so its other
  alignments in the other candidate are kept as hidden rows (`region.LAYER_PARTIAL`,
  `_partial_alignments`, `_partial_lines`): those that do not qualify at the
  display pid and plen but reach "partial pid" (default 0.99) with plen from
  "partial plen" (default 0.3) -- partial ones of any type (a clipped or split
  primary, a supplementary, a secondary), and whole-length ones short of the
  display pid that hold an indel of at least "min indel" (one without is
  divergence, not an SV) -- inside the locus, each a line from the read's whole
  placement. Only a read's best of them are kept, as for placements: those within
  best_pid_delta of its best pid there (at most 5), with the supplementary pieces
  of a kept primary; without that, a read in a repeat brought every copy
  (atg003038l_translocation_002: 471,930 rows and a 575 s build at partial pid
  0.99; 13,076 rows and 46 s now, against 44 s without them).
  They are never drawn unselected; when a marginal-bin hover or a tap selects the
  read they show in orange (`sel_alpha`). The tooltip says what the alignment is,
  its pid and plen, whether a clip runs off a contig end (within 1 kb,
  `CLIP_END_TOL`) or lies inside a contig, and its indels of at least "min indel"
  (default 50 bp), each also marked with a cross on the line, shown with it
  (`RegionData.indels`, `glyphs._draw_indels`). The read keeps its class, the
  marginals do not count these alignments, and the "partial" layer checkbox
  switches them and their crosses off.

  Why a separate pid and the indels: of the one-sided reads that had no partial
  row at the display pid, most had alignments in the other candidate short of it.
  A piece ending at a breakpoint runs past it before the aligner clips, which
  costs a few tenths of a percent of pid (ptg000114l_translocation_001: median
  0.988); and an SV of a few hundred bp is often kept inside one whole-length
  alignment as a single long gap (atg003038l_translocation_002: 23% of such
  alignments hold an indel of at least 50 bp, 9% one of 500 bp). Pieces that are
  divergent throughout -- a paralog hit where the read's own sequence is missing,
  median pid 0.945 on atg004724l_translocation_001 -- stay below 0.99.
- **What a control costs.** Only pid, plen, partial plen, partial pid, min indel, best_pid_delta, pid_tie_eps and the
  bin size change which reads are drawn; they, and a resize, ask the server for a
  new view. Colouring, C, layers, the reads mode and the background toggle are
  applied by the page to the plot it holds (`restyle` in `static/app.js`): the
  line rows carry their layer, class colour and read, the read table
  (`qv-reads`) each read's log-ratio, and both marginal stackings are drawn with
  one shown (renderers tagged `stack:class` / `stack:ratio`, the hover reading the
  hidden `qv-mode` Div). `lineStyle` mirrors `glyphs.line_style`, so a restyle
  draws exactly what the server would (pixel-identical in the browser check), and
  C also reworks the score card from the per-read log-ratios the view carries.
  Every per-line column and callback index is a numpy array, which Bokeh sends as
  one binary buffer: names, classes and contigs are integer codes resolved by the
  tooltip (`CustomJSHover`), colours packed uint32. On a 70,000-line locus a view
  of an already-built region went from about 6 s and 12 MB to 0.6 s and 2.1 MB.
  Building a region for new thresholds (about 5 s on 6,000 reads) is what is left.
- **Zooming one axis.** The wheel over the dotplot zooms both axes; over a
  marginal it zooms that marginal's genome axis alone -- the top one X, the right
  one Y -- which the dotplot shares, so the dotplot follows on that axis only. The
  marginals' coverage axis, one scale for both, is never zoomed by the wheel.
- **Marginals on one scale.** As in `report_plot`, the two coverage marginals
  are the same size and share one coverage range, 0 to the larger of their two
  peaks (one `Range1d`, so zooming one zooms the other): equal coverage is
  equally long on either, and the candidates compare at a glance. The peak is
  of what is in view, so the scale follows the zoom.
- **Marginals rebinned on zoom.** The view carries each read's best placements
  as intervals on each axis (`RegionData.placements`), and the page bins them
  afresh for the range in view on every zoom or pan, at most once a frame
  (`glyphs.REBIN_JS`): about 1000 bins across the view, the bin size rounded up
  to 1, 2 or 5 x 10^k bp, never below 1 bp, bins restarting at each segment.
  A placement adds its weight to every bin it touches, by a difference array per
  band, so a rebin is one pass over the placements -- about 11 ms for both
  marginals on 200,000 placements; the dotplot's redraw costs far more. As bins
  shrink a bar tends from "reads touching the bin" to depth. A hovered bin finds
  its reads by scanning the placements, a selected read its bins by binary search,
  so neither needs a stored bin index. `region.marginal_bins` is the Python twin,
  drawing the first frame; `tests/test_rebin.py` runs both and checks them against
  the server's fixed-size bins. The `binsize` control now sets only the stretch
  a shift+click aligns.
- **Sorting along the diagonal.** The "sort" control lays one axis's segments
  along the other's diagonal (`region.diagonal_order`): a segment goes where the
  weighted median of its reads' lines lands on the other axis, each line weighing
  its length over its read's line count (a repeat read pulls no harder than a
  unique one; partial rows do not vote), and it is flipped when those lines
  mostly run backwards. Segments no line reaches go last. `region.arrange` then
  moves every axis coordinate with its segment -- lines, placements, indel marks
  -- with no new fetch; the region is cached per sort. A flipped segment runs
  from its end to its start: everything mapping the axis to the genome reads the
  flags (`RegionData.flips`) -- the readout (marked "rev"), marginal bins, contig
  end marks, the rectangle aligner, remembered rectangles (placed after the view
  loads, in its layout) -- and the k-mer background is resampled into the new
  layout to within one of its bins. A segment whose reads land in two places on
  the other axis (a duplication) goes by the heavier one.
- **Coverage dropouts.** Each marginal marks in thin red every stretch of a
  segment that no read's best placement covers (`region.coverage_gaps`), exact
  to the bp and drawn at least a pixel wide, so a dropout narrower than a bin
  still shows when zoomed out.
- **Resizing by the borders.** Four handles on panel borders resize panels in
  place -- sizes set on the page (`applySizes`, finding panels by their `qv:`
  tags), no new view: the dotplot's right border sets its width and its bottom
  border its height (independently, so its pixels need no longer be square); the
  border between the top marginal and the dotplot, and the right marginal's outer
  edge, set the marginals' size -- both together, so their shared coverage scale
  stays comparable. A drag shows the border moving and applies on release, so a
  large locus is not redrawn on every mouse move. The contig-end triangles read
  their bar's current length. Sizes chosen this way are kept across views and
  loci until a handle is double-clicked, which goes back to the sizes fitted to
  the window; while they are kept, resizing the window does not refit the plot.
- **One plot at a time.** A view embeds into a container of its own and is
  numbered; one that finishes after a newer one was asked for removes itself, and
  the newest leaves no other Bokeh document. Before, only the fetch was aborted: a
  view still embedding when a newer one arrived stayed beside it, and two stacked
  plots added their translucent lines up -- the plot brightened after changing
  thresholds.
- **Contig ends.** On the contig bars beside the dotplot, a black triangle marks a
  segment edge that is the start or the end of its contig, pointing into the
  segment. An edge within 1 kb of the contig's end counts (`CONTIG_END_TOL`):
  reads seldom align to a contig's last few hundred bases, so recruited segments
  stop just short; the hover text gives the contig's length and how far past the
  segment its end lies. Lengths come from the BAM headers
  (`RegionData.contig_lengths`). The triangles keep their size in pixels and sit
  half of it inside their segment at any zoom. `glyphs.contig_end_marks`.
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
| `GET /api/events` | server-sent events: queue progress and row updates (a locus re-recruited away comes as `{id, removed, moved_to}`) |
| `GET /api/loci/{id}/seeds` | the locus's seeds with the recruitment parameters each was recruited with |
| `POST /api/recruit` | `{seeds, params}`: recruit those seeds again with those parameters |
| `GET /api/loci/{id}/view?w=&h=&...` | Bokeh `json_item` of the report layout for that size and those settings, plus the score summary |

Launch: `venv/bin/python -m qaview [project.toml] [--port 5006]`. Without a
project file the page starts on a form in the left sidebar: the project file's
path on the server (a browser upload would break the paths inside it, which are
relative to the file), with the last thirty projects opened one click away
(`~/.config/qaview/recent.json`). "change" beside the project name brings the
form back; opening another project stops the current one's recruitment and
scoring and reloads the page onto the new one (`POST /api/open`; everything
that belongs to a project lives in a `Session` replaced whole).

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
