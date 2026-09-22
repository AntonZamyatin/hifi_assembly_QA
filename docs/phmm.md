# Pair HMM for read mark likelihood — development plan

Implements §5 of `model.tex`: the mark likelihood `e(Y | p, A)`, replacing the
`exp(-NM)` placeholder in `plotters/likelihood.py`.

Scope is the emission term only. The count/exposure model (`-2λ L_eff`), the cap
`C`, and the tolerance `η` are separate work.

Sizing is relative, not calendar: **S** mechanical, **M** real work but the shape
is known, **L** can go sideways, **XL** a project of its own.

## Design decisions taken up front

**One core function, parameterised by window.** The primitive is

```python
forward_window(read, qual, template, params, origin=None) -> WindowScore
```

returning `log Σ_{p ∈ origin} P(r, ℓ | Q, p)` plus diagnostics. A thin adapter
extracts the template from a FASTA given `(contig, pos, pad)` — no database, no
recruiter, no BAM. Simulation and real data hit the same code path.

**A window has two extents, and they are not the same interval.** Its *template*
extent must be long enough to hold the read plus indel slack. Its *origin*
extent is the set of start positions the window is responsible for. Conflating
them is a correctness bug, so `origin` controls the second independently:
`None` for the whole template, `(lo, hi)` for a half-open range, an integer for
one position.

This matters because the model's `Σ_p` runs over every origin in the recruited
*unit*, while a padded window around one alignment covers only that alignment's
neighbourhood. A read in a repeat has several placements — and those secondary
placements *are* the copy-number signal — so in general one call is one term of
the sum. Terms combine by log-sum-exp, valid only when the origin ranges are
disjoint: template extents overlap by design, origin extents must not, or shared
origins are counted twice and every repeat is silently inflated. The whole-unit
free-start pass (`origin=None` over the unit) is the case where one call is the
entire sum, and it is the oracle stage F is measured against.

**Padding is a parameter with a self-check.** The true origin differs from a
reported alignment start by indel drift plus soft-clipping, so the window must
be `ℓ + 2·pad` with pad covering that offset. Too small truncates the origin sum
silently. `WindowScore` therefore reports the fraction of terminal mass reaching
the window edges; a score is only trusted when edge mass is negligible *and*
stable under doubling pad. This makes the API self-validating, which matters
because every later stage depends on it.

Pad should be per read, not a global constant. The `aln` table already carries
`soft_left`, `soft_right`, `hard_left`, `hard_right`, and a read reported at
`pos` with `soft_left = c` has its true origin near `pos − c`. Sizing from those
plus slack for indel drift is cheaper than a uniform large pad, and the edge-mass
diagnostic catches the cases where it was still too small.

**Context-independent gaps first.** With `(δ_I, δ_D, γ_I, γ_D, τ)` constant, the
normaliser `Z(ℓ, p, A)` depends on neither `p` nor `A` (§5.6), so it cancels from
ΔLLR — and `τ` cancels with it. The two-tier scheme of §5.8 collapses to a single
pass. Context dependence is deferred and must be justified by measurement.

**`τ = 0`: read length is conditioned on, not generated.** A HiFi read ends
because the size-selected SMRTbell insert ran out, not because a per-base
stopping process fired, so the geometric length distribution a positive `τ`
implies (mean `1/τ`) describes no real library — as §5.6 concedes. The model
discards it anyway by dividing through `Z(ℓ)` and letting the empirical `f(ℓ)`
carry length, which makes `τ` a normalisation device rather than a parameter.
It must never be fitted to observed read lengths.

`forward_window` therefore returns `log Σ_m [f_M(ℓ,m) + f_I(ℓ,m)]` and leaves
the termination prefactor off. That prefactor is uniform — every path with `ℓ`
emitting states terminates exactly once, from M or I, at probability `τ` either
way — which is why §5.5 eq. (9) carries it out front. What does not factor out
is `τ` inside `M→M` and `I→M`, since paths with different indel counts carry
different powers of those; measured at ℓ = 15 kb that is a −1.5 nat shift on the
read but only −0.01 nats between candidates disagreeing by 100 indels, so it is
common mode. A positive `τ` remains available for the joint model over `(r, ℓ)`.

**`B→D` and `D→E` stay excluded, and padding is the reason.** A padded window
invites adding them to walk across the flanks. They are not needed: the free
start absorbs leading padding (the first read base may align at any template
position) and the free end absorbs trailing padding (the sum runs over every
`m`). Adding them on top would count "begin at `v`, delete `t_v`" and "begin at
`v+1`" as two paths for one alignment, and symmetrically at the end. A classical
pair HMM needs `B→D` because it starts at a fixed position; a free-start pass has
already spent that degree of freedom on the origin sum, which is where the
copy-number signal lives. Remark 3 attributes the degeneracy correctly — it comes
from marginalising over origins — so a wider window makes the exclusions more
necessary, not less.

**The full-table engine is permanent.** `_dp_full_table` in `dp.py` is the oracle
every faster path is tested against, not scaffolding. It keeps the scalar
deletion loop for exactly that reason.

## Module layout

```
src/qamodel/
  exposure.py      lambda from genome-wide depth; L_eff; exposure term   [I, done]
  score.py         per-read ratio, clamp at +/-C, Delta S; score_locus   [I, done]
  phmm/
    params.py      PhmmParams; transitions; quality table                [A, done]
    dp.py          the recursion: full-table oracle + two-row engine     [B, F1]
    band.py        which cells to compute; diagonal or CIGAR-anchored    [F2]
    kernel.py      compiled engine, linear space with rescaling          [F4]
    sample.py      the model run forwards                                [D]
    reads.py       read sequence and qualities from primary records      [H]
    windows.py     placements -> strand-oriented windows, disjoint origins [F3]
```

`qamodel` imports nothing from `qatool`; the dependency runs one way only
(`qatool.plotters.likelihood → qamodel`), which is what keeps stages A–E usable
on simulated sequence with no database, recruiter or BAM in play.

## Definitions

Stated strictly, because the difference between them is where the errors are.

**Unit** `U_A` — a finite set of disjoint intervals of assembly `A` (the
recruited segments), possibly on different contigs. Total length `L_A`.

**Origin** `p = (c, s, o)` — contig, 0-based position, strand `o ∈ {+,−}`, with
`s` inside a segment of `U_A`. The slot set has exactly `2·L_A` members, which is
the `2L_A` of §7's uniform prior and of the exposure `−2λL_A`.

**Template at `p`** — a function of `p` alone, with `[b, e)` the segment holding
`s`: `A[c][s:e]` when `o = +`, `revcomp(A[c][b:s+1])` when `o = −`. Truncated at
the segment boundary, so origins near an edge get a short template — which *is*
Proposition 8's containment problem, not a separate concern.

**Mark** `Y = (r, Q)` — the read as sequenced, qualities in sequencing order,
never reverse-complemented.

**`e(Y | p, A)`** — one `forward(r, Q, template_at_p, origin=0)` call.

An alignment is **not** part of any of these. Alignments are a device for finding
origins with non-negligible `e`; an origin the aligner never reported is still a
term in `Σ_p`, and since Proposition 4's signal is the *count* of origins,
truncating the sum at `-N` attacks the quantity being measured (§12.4).

## Stages

Stages A–E need no project data at all. Real-data plumbing is stage H.

### Stage A — Parameters and scalars (S) — DONE

`PhmmParams` dataclass, defaults, quality-to-`ε` conversion, `s = log(3(1-ε)/ε)`.

Calibration anchor: this data is mean Q39, so ε ≈ 1.2×10⁻⁴ and `s ≈ 10.1` nats.
The current `exp(-NM)` uses `s = 1`.

*Accepts when:* `s` matches the doc's table (5.7 @ Q20, 8.0 @ Q30, 10.3 @ Q40).

### Stage B — `forward_window` and Viterbi (M) — DONE

numpy, log space with LSE, no banding, small inputs. Implements §5.5 exactly,
including the two forbidden transitions (`B→D`, `D→E`).

**Strand lives in the template, not the read.** §5.2 defines `t` as the assembly
sequence read from `s` on strand `o`, so a reverse-strand origin takes
`revcomp(A[c][b : s+1])` as its template while the read stays as sequenced with
its qualities in sequencing order. Reverse-complementing the read instead is
*not* equivalent: measured, it agrees to 1e-6 on an ordinary read but differs by
0.25 nats when the read begins with an insertion, because `B→I = δ_I` and
`B→M = 1−δ_I` are asymmetric and reversing a path swaps which end pays which.
Nothing in `dp.py` reverse-complements anything; the caller supplies the
strand-oriented template, and stage H's adapter is where that happens.

*Accepts when:*
- **Exact path enumeration.** For ℓ ≤ 4, m ≤ 6, brute-force every state path and
  match to 1e-12. Ground truth, and the most valuable test in the suite.
- **Window additivity.** `forward_window` over a template equals `logsumexp` of
  pad-0 calls at each origin, to 1e-9. This is the definition of correctness for
  the windowing.
- **Windows compose over disjoint origin ranges.** Partitioning a unit's origins
  across several calls, each with a template running to the unit's end, must sum
  to the whole-unit pass to 1e-12 — the rule stage F's anchored kernel depends on.
  Overlapping the ranges instead must visibly inflate the total, since that is
  the silent failure mode it guards against.
- **Homopolymer degeneracy.** Template `AAAAAAA`, read `AAAAAA`: forward/Viterbi
  = exactly **5**, not 7. A deletion at template position 1 needs `B→D` and at
  position 7 needs `D→E`, both forbidden; the 5 survivors carry an identical
  transition product (4 `M→M` each, independent of position). model.tex Remark 6
  says 7, which is the count before those exclusions — the doc needs correcting.
  The test is better at 5 because it exercises both.
- **Mismatch cost.** Each substitution costs `s` within 1e-9.
- **Copy number.** A template with 3 exact tandem copies scores `log 3` above one
  copy, within 1e-6. Proposition 4's mechanism in miniature; standing regression
  test, because if it breaks the method is blind to collapses.
- **Edge mass.** Reported edge mass falls monotonically as pad grows, and the
  score converges.
- **Termination is a prefactor.** At `τ = 0` the score is finite and the
  termination edges are never consulted; at `τ > 0` the score shifts by exactly
  `ℓ`-independent `log τ` relative to `τ = 0`, confirming eq. (9) carries it out
  front and that no path terminates twice.
- **No origin double-counting.** With gaps disabled, padding a template by `k`
  bases on either side must leave the score unchanged to machine precision: each
  origin is counted once, and the added origins explain a random read only
  through ~L mismatches. With gaps enabled the score *does* rise, and should —
  a pad that extends a homopolymer across the boundary lets the read's last base
  match further right after deletions, worth `(1 − γ_D^k)·δ_D`, which is
  Remark 6's degeneracy in a second guise. The increment is data-dependent but
  non-negative and bounded by the total gap-open mass
  `δ_D/(1−γ_D) + δ_I/(1−γ_I)`.

### Stage C — Normaliser Z (S)

Same recursion, all emissions set to 1.

*Accepts when:*
- **Exact.** For ℓ = 4 enumerate all 4⁴ = 256 reads, sum, match `Z(4)` to 1e-12.
  Emissions sum to one per state, so this holds exactly or the code is wrong.
- **Template independence.** `Z(ℓ)` identical across unrelated templates —
  confirming it cancels and can be skipped in production until context
  dependence lands.

### Stage D — Generative sampler (S) — DONE

`sample_read(template, params, rng) -> (read, qual, true_path, true_origin)`.
Walk the chain, emit, record the path. Small, and it unlocks every stage below.

*Accepts when:*
- **Empirical/analytic agreement. Met.** With `τ > 0` the chain terminates on its
  own and a read's frequency is exactly `τ·exp(forward)`. Over 200k draws the
  worst of the top 20 reads is 2.4 standard errors out, and summed over all 84
  reads of length ≤ 3 the observed mass is 0.5787 against 0.5781 predicted.
- **Length-conditioned agreement. Met.** With `τ = 0` and a fixed length the
  sampler draws from `forward / Z(ℓ)`; normalising both sides removes the unknown
  constant, and the worst gap is 0.0012. Inverting it estimates `Z(6) ≈ 1.003`,
  which stage C should reproduce directly.
- **Indel lengths geometric. Met.** Insertion runs mean 2.228 against 2.222;
  deletion runs 1.538 against 1.538.

Two test premises were wrong and are recorded because both recur. Selecting reads
by observed count biases a global mass check upward — it keeps exactly the reads
that came up more often than expected — so the check enumerates every read of
length ≤ 3 independently of the sample instead. And discarding truncated walks
censors deletion runs, because a D run still going when the template runs out is
what causes truncation; on a 40 bp template that pulled the deletion mean 4
standard errors low while leaving insertions untouched, since they consume no
template. The test now uses a 400 bp template, where nothing truncates.

### Stage E — Propositions on simulated data (M)

The point of this stage: check the doc's closed forms against the implementation
before real data can confuse the picture. Each is a synthetic pair of templates,
reads sampled from the truth, scored under both.

| construct | prediction | source |
|---|---|---|
| tandem array, k₁ vs k₂ copies | per-read gain `log(k₁/k₂)` | Prop. 4 |
| false duplication | per-read `−log 2` | Prop. 6 |
| two haplotypes at divergence d | per-read `log((1+u)/u) ≈ d·ℓ·s` | §10 |
| collapse junction reads | one-sided, ratio ≫ 0 | Remark 12 |

*Accepts when:* each measured mean matches the closed form within Monte Carlo
error. Disagreement is informative in both directions — it finds bugs in the code
*and* errors in the doc, as the `log 5` case already shows.

Realistic simulation (pbsim3 HiFi profiles) belongs after this, to probe
misspecification. Sampling from the model itself comes first because only there
is the truth exact.

### Stage F — Fast engines and origin enumeration (L) — F1, F2, F4 DONE

The risky stage. Origin candidates from a per-locus minimap2 run against the
extracted unit (`-N 200 -p 0.05`). Banded DP anchored on each candidate's
diagonal, adaptive width, band-edge mass reported as a diagnostic (§5.8).

Why anchored rather than free-start over the whole unit: free-start on a 100 kb
unit with a 15 kb read is ~1.5×10⁹ cells, ×3 states ×~500 reads ×2 candidates
≈ 4×10¹² cell updates per locus — fine for curated loci, hopeless for a survey.
Anchoring costs ~6×10⁷ cells per read at k=20, w=200.

§12.4 objects that aligner secondary limits truncate the origin sum and kill the
copy-number signal. That only bites genome-wide: re-aligning recruited reads
against the extracted ~100 kb unit with generous settings costs seconds and
sidesteps the global `-N 20` entirely.

**F1 (done)** row recycling and a vectorised deletion row, in `dp.py`. Memory
went from 5.1 GB to 9 KB for a 15 kb read, which is what made real reads possible
at all.

**F2 (done)** banding, in `band.py`. Diagonal or CIGAR-anchored, with `band_edge`
reporting mass pressed against a band boundary. 8.8× over unbanded.

**F4 (done)** the compiled engine, in `kernel.py`. **28.5× over banded numpy**,
so about 13 s per locus for 500 reads against two assemblies rather than 6 min.

The linear-space question was settled by measurement, and the answer reverses
what this document said before. At band width 200 `np.logaddexp` costs 3.66 µs
per call against 0.41 µs for `np.maximum`, and ~89% of that is the `exp`/`log1p`
pair; across a row those two operations are ~85% of the time. In linear space
with per-row rescaling they vanish — every cell is multiply-add, ~13 flops. That
is where the 28× comes from, not from compilation as such. Log space remains
right for the numpy engine, where the deletion scan must go through
`logaddexp.accumulate` because a Python loop is unaffordable; inside numba a
sequential loop is exactly what is wanted.

**F3 (done)** `phmm/windows.py`. Containment on the aligned span, the recruiter's
own rule; a read crossing a segment edge in either candidate is dropped from both.
Clips are in reference orientation, so a reverse-strand read starts at
`end + clip_right - 1`, with the template reverse-complemented and the CIGAR
reversed. Placements whose origin ranges would overlap merge into one window
owning the hull; the band is the union of members' CIGAR bands. Checked on
synthetic sequence for both strands and both clip sides, and three tandem copies
give `log 3` through the full placement-to-window path.

*Accepts when:*
- Matches `dp.py` to 1e-9 on every stage B–D case. **Met**: 142 configurations
  for the kernel, worst 8.5e-14; terminal rows agree to 6.3e-13.
- **The budget experiment.** On a real repeat locus, anchored `Σ_p` matches the
  free-start oracle within 0.1 nats. Half met — disjoint windows already
  reproduce the whole-unit pass to 2e-16 on synthetic sequence; the real-locus
  half belongs to F3, with real aligner placements.
- Band-edge mass below a stated threshold; throughput recorded. **Met.**

### Stage G — Parameter estimation (M)

Method of moments from CIGARs of high-identity primary alignments in undisputed
sequence: gap opens and extends give `(δ_I, δ_D, γ_I, γ_D)`. Fitted once, held
fixed across both candidates (§12.5). `τ` cancels under context-independence.

*Accepts when:*
- **Round trip.** Parameters re-estimated from reads sampled by stage D recover
  the generating values within Monte Carlo error. This is an identifiability
  check and is only possible because the sampler exists.
- No per-candidate fitting path exists in the API — the constraint is structural,
  not a convention.

Baum–Welch is the better estimator but waits until context dependence makes `Z`
matter.

### Stage H — Real-data adapters (M) — DONE

`ReadStore` and locus-level origin enumeration.

Prerequisite discovered while planning: secondary alignment records carry no
`SEQ`/`QUAL` (SAM spec, verified on `hg002.hifiasm_hic.sorted.bam`), and at a
repeat locus the multiplicity reads appear mostly as secondaries. No raw HG002
HiFi FASTQ is on disk. Verified route: the `aln` table indexes every placement by
`(assembly, read_name)` with `aln_type`, so look up the read's `primary` row,
fetch that locus in the BAM, take `query_sequence` / `query_qualities`. A locus
needs ~500 such lookups — seconds, and cacheable.

*Accepts when:* returned length matches `q_len` in the DB across all four
`aln_type` values; a read whose only local placement is secondary still resolves;
a missing primary raises rather than returning silently.

### Stage I — The full locus score (M) — DONE

Built as the whole score rather than the emission term alone, because an
emission-only number prefers a false duplication by exactly `log 2` per read
(Proposition 6, reproduced in `tests/test_score.py`).

    Delta S = sum_j clamp(log sum_A1 e_jp - log sum_A2 e_jp, -C, +C) - 2 lambda (L_eff1 - L_eff2)

`qamodel/exposure.py`: `lambda = N / (2 L_dip)` from primary reads genome-wide —
1.05e-3 on HG002, 32x per haplotype — and `L_eff` in the exact form over the
empirical read-length histogram. `qamodel/score.py`: every read clamped at `±C`
(default 200), so the score is continuous and a one-sided read simply lands on
the clamp. `LocusScore.at(C)` re-evaluates for another C from stored ratios.
Entry point `qatool.plotters.phmm_likelihood_report`.

**HG002 seed locus:** Delta S = −271.95 nats for `curated` over `hifiasm_hic`,
106 reads, none dropped, 28 s. The two candidates are identical here except for
one event — an (AT)n microsatellite three units shorter in `hifiasm_hic`, against
an Alu poly-T tail at `h2tg000020l:61,904,626`. All 18 reads that favour
`curated` span it; no read favours `hifiasm_hic` by more than 0.69 nats.
Regression test: `tests/test_score_hg002.py`.

**Caveat on magnitude.** The event is a microsatellite length difference — the
context where HiFi indel error concentrates and where context-independent gap
parameters are least right. The *direction* rests on 18 spanning reads against
none; the *size* depends on uncalibrated placeholders. Stage G, and then the
deferred context-dependent gaps, are what make the magnitude interpretable.

## Deferred

- **Context-dependent gaps** (§5.7) — **XL**. ~40 parameters, and it resurrects
  `Z` and the two-tier scheme. Justify by measuring whether homopolymer
  misspecification moves verdicts; stage E is where that gets measured.
- **Per-origin normalisation** (tier two) — **M**, and only meaningful once
  context dependence lands.
