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
forward_window(read, qual, template, params) -> WindowScore
```

returning `log Σ_p P(r, ℓ | Q, p)` over every origin `p` in `template`, plus
diagnostics. Everything is a special case: pad = 0 is fixed-origin, a 100 kb
template is the free-start pass of model.tex §5.8. A thin adapter extracts the
window from a FASTA given `(contig, pos, pad)` — no database, no recruiter, no
BAM. Simulation and real data hit the same code path.

**Padding is a parameter with a self-check.** The true origin differs from a
reported alignment start by indel drift plus soft-clipping, so the window must
be `ℓ + 2·pad` with pad covering that offset. Too small truncates the origin sum
silently. `WindowScore` therefore reports the fraction of terminal mass reaching
the window edges; a score is only trusted when edge mass is negligible *and*
stable under doubling pad. This makes the API self-validating, which matters
because every later stage depends on it.

**Context-independent gaps first.** With `(δ_I, δ_D, γ_I, γ_D, τ)` constant, the
normaliser `Z(ℓ, p, A)` depends on neither `p` nor `A` (§5.6), so it cancels from
ΔLLR — and `τ` cancels with it. The two-tier scheme of §5.8 collapses to a single
pass. Context dependence is deferred and must be justified by measurement.

**The reference implementation is permanent.** `reference.py` is the oracle every
later stage is tested against, not scaffolding.

## Module layout

```
src/qatool/phmm/
  params.py      PhmmParams; defaults; estimation from CIGARs
  reference.py   numpy oracle: forward_window, viterbi_window, normaliser
  sample.py      generative sampler: draw reads from the model
  window.py      FASTA adapter: (contig, pos, pad) -> template
  kernel.py      numba banded kernels
  origins.py     candidate-origin enumeration per locus
  reads.py       ReadStore: read_name -> (seq, qual)
  score.py       per-read log Σ_p e_jp and the log-ratio
```

## Stages

Stages A–E need no project data at all. Real-data plumbing is stage H.

### Stage A — Parameters and scalars (S)

`PhmmParams` dataclass, defaults, quality-to-`ε` conversion, `s = log(3(1-ε)/ε)`.

Calibration anchor: this data is mean Q39, so ε ≈ 1.2×10⁻⁴ and `s ≈ 10.1` nats.
The current `exp(-NM)` uses `s = 1`.

*Accepts when:* `s` matches the doc's table (5.7 @ Q20, 8.0 @ Q30, 10.3 @ Q40).

### Stage B — `forward_window` and Viterbi (M)

numpy, log space with LSE, no banding, small inputs. Implements §5.5 exactly,
including the two forbidden transitions (`B→D`, `D→E`). Reverse-strand reads
handled by reverse-complementing the read, once, at the adapter boundary.

*Accepts when:*
- **Exact path enumeration.** For ℓ ≤ 4, m ≤ 6, brute-force every state path and
  match to 1e-12. Ground truth, and the most valuable test in the suite.
- **Window additivity.** `forward_window` over a template equals `logsumexp` of
  pad-0 calls at each origin, to 1e-9. This is the definition of correctness for
  the windowing.
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

### Stage C — Normaliser Z (S)

Same recursion, all emissions set to 1.

*Accepts when:*
- **Exact.** For ℓ = 4 enumerate all 4⁴ = 256 reads, sum, match `Z(4)` to 1e-12.
  Emissions sum to one per state, so this holds exactly or the code is wrong.
- **Template independence.** `Z(ℓ)` identical across unrelated templates —
  confirming it cancels and can be skipped in production until context
  dependence lands.

### Stage D — Generative sampler (S)

`sample_read(template, params, rng) -> (read, qual, true_path, true_origin)`.
Walk the chain, emit, record the path. Small, and it unlocks every stage below.

*Accepts when:*
- **Empirical/analytic agreement.** Over many samples from a short template, the
  observed frequency of each read matches `exp(forward_window)` — the sampler and
  the scorer are then provably the same model.
- Observed indel-length distribution is geometric with the specified `γ`.

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

### Stage F — Banded numba kernel and origin enumeration (L)

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

Consider row-rescaled linear space over LSE — avoids an `exp` and a `log` per
cell, typically 3–5× faster — but only if the stage B/C tests still pass to 1e-9.

*Accepts when:*
- Matches `reference.py` to 1e-9 on every stage B–D case.
- **The budget experiment.** On a real repeat locus, anchored `Σ_p` matches the
  free-start oracle within 0.1 nats. This gate sets the compute budget for
  everything downstream and must be settled here, not after wiring in.
- Band-edge mass below a stated threshold; throughput recorded in cells/s.

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

### Stage H — Real-data adapters (M)

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

### Stage I — Wire into the score (S)

Replace `exp(-NM)` in `plotters/likelihood.py` with `log Σ_p e_jp`. Emission term
only — leave `-n log L_A` alone so the two changes stay attributable.

*Accepts when:* on the HG002 seed locus, per-read log-ratios are finite, read
ranking is stable under a 2× change in band width, and wall time is recorded.

## Deferred

- **Context-dependent gaps** (§5.7) — **XL**. ~40 parameters, and it resurrects
  `Z` and the two-tier scheme. Justify by measuring whether homopolymer
  misspecification moves verdicts; stage E is where that gets measured.
- **Per-origin normalisation** (tier two) — **M**, and only meaningful once
  context dependence lands.
