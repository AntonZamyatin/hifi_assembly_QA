"""Stage G: gap parameters from aligned reads (docs/phmm.md).

The four gap parameters are the transition probabilities of the pair HMM, so
the estimate is a transition count along an alignment path:

    delta_I = n(M->I) / n(M->.)      gamma_I = n(I->I) / n(I->.)
    delta_D = n(M->D) / n(M->.)      gamma_D = n(D->D) / n(D->.)

That is the maximum-likelihood estimate when the path is known, and the path a
CIGAR gives is the aligner's, not ours.  Two consequences are worth stating
rather than burying:

**The aligner's path is not the model's posterior.**  minimap2 chooses gaps
under its own affine penalties (`map-hifi`: O 4,24 / E 2,1), so a gap it would
not pay for is invisible here, and where several placements of one gap tie it
picks one convention.  The estimate is therefore of "gaps as minimap2 represents
them".  `viterbi_counts` re-counts the same reads under our own model so the two
can be compared; the round trip on sampled reads (`tests/test_estimate.py`)
checks the estimator itself, where the path really is the model's.

**Long gaps are not this model's gaps.**  A 138 bp deletion in a HiFi read is a
variant, an assembly error, or a collapsed repeat -- not a sequencing indel, and
including it inflates gamma_D towards 1.  Runs longer than `max_gap` split the
path: the long run and the two transitions touching it are dropped, and the
alignment continues as two independent paths.  Dropping the whole read instead
would bias the other way, since it is the reads crossing hard sequence that
carry most short gaps too.  `Counts.split` records how often this fired.

The transitions I->D and D->I do not exist in the model (params.py); CIGARs do
contain them, so they are counted separately and reported as
`Counts.unmodelled` rather than being folded into a neighbour.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass, field, asdict

import numpy as np

from qamodel.phmm.params import PhmmParams, eps_from_phred

MAX_RUN_HIST = 20           # run lengths tabulated individually; longer ones pool
MAX_PHRED_BINS = 94         # SAM caps Phred at 93
DEFAULT_MAX_GAP = 20        # a run longer than this is not a sequencing indel


@dataclass
class Counts:
    """Transition and base counts over many alignment paths."""

    m_m: int = 0
    m_i: int = 0
    m_d: int = 0
    i_i: int = 0
    i_m: int = 0
    d_d: int = 0
    d_m: int = 0
    i_d: int = 0                # unmodelled adjacency, counted not folded
    d_i: int = 0
    m_bases: int = 0
    i_bases: int = 0
    d_bases: int = 0
    mismatches: int = 0         # from NM, where available
    eps_sum: float = 0.0        # sum of per-base eps over aligned read bases
    eps_bases: int = 0
    reads: int = 0
    split: int = 0              # paths split at an over-long run
    # Attrition, so the kept set can be seen rather than assumed.  Every filter
    # selects on alignment quality, which is correlated with the thing being
    # measured: keeping too little biases the gap rate down.
    skip_mapq: int = 0
    skip_split_read: int = 0
    skip_clip: int = 0
    skip_identity: int = 0
    i_hist: list = field(default_factory=lambda: [0] * (MAX_RUN_HIST + 1))
    d_hist: list = field(default_factory=lambda: [0] * (MAX_RUN_HIST + 1))

    @property
    def unmodelled(self) -> int:
        return self.i_d + self.d_i

    def __iadd__(self, other: "Counts") -> "Counts":
        for k, v in vars(other).items():
            if isinstance(v, list):
                mine = getattr(self, k)
                for i, n in enumerate(v):
                    mine[i] += n
            else:
                setattr(self, k, getattr(self, k) + v)
        return self


# pysam CIGAR op codes.  M/=/X all consume both and are one state here.
_MATCH, _INS, _DEL, _REF_SKIP, _SOFT, _HARD, _PAD, _EQ, _DIFF = range(9)
_STATE = {_MATCH: "M", _EQ: "M", _DIFF: "M", _INS: "I", _DEL: "D"}


def count_path(cigartuples, counts: Counts, max_gap: int = DEFAULT_MAX_GAP) -> None:
    """Add one alignment's transitions to `counts`.

    `cigartuples` is pysam's (op, length) list.  Clips and reference skips end
    the current path rather than transitioning through it: a clip means the read
    left the template, which the free start and free end already handle.
    """
    prev = None                                  # state of the path's last base
    for op, n in cigartuples:
        state = _STATE.get(op)
        if state is None or n == 0:              # clip, skip, pad: break the path
            prev = None
            continue
        if state != "M" and n > max_gap:         # variant-sized: not our gap model
            counts.split += 1
            prev = None
            continue

        if state == "M":
            counts.m_bases += n
            counts.m_m += n - 1
        elif state == "I":
            counts.i_bases += n
            counts.i_i += n - 1
            counts.i_hist[min(n, MAX_RUN_HIST)] += 1
        else:
            counts.d_bases += n
            counts.d_d += n - 1
            counts.d_hist[min(n, MAX_RUN_HIST)] += 1

        if prev is not None:
            setattr(counts, f"{prev.lower()}_{state.lower()}",
                    getattr(counts, f"{prev.lower()}_{state.lower()}") + 1)
        prev = state


@dataclass(frozen=True)
class Estimate:
    """Estimated parameters with their binomial standard errors."""

    params: PhmmParams
    se_delta_i: float
    se_delta_d: float
    se_gamma_i: float
    se_gamma_d: float
    counts: Counts

    def summary(self) -> str:
        c, p = self.counts, self.params
        gaps = c.i_bases + c.d_bases
        mm_rate = c.mismatches / c.m_bases if c.m_bases else float("nan")
        eps = c.eps_sum / c.eps_bases if c.eps_bases else float("nan")
        return (
            f"{c.reads:,} reads, {c.m_bases:,} aligned bases\n"
            f"  delta_i {p.delta_i:.3e} +/- {self.se_delta_i:.1e}   "
            f"gamma_i {p.gamma_i:.4f} +/- {self.se_gamma_i:.4f}  "
            f"(mean run {1 / (1 - p.gamma_i):.2f})\n"
            f"  delta_d {p.delta_d:.3e} +/- {self.se_delta_d:.1e}   "
            f"gamma_d {p.gamma_d:.4f} +/- {self.se_gamma_d:.4f}  "
            f"(mean run {1 / (1 - p.gamma_d):.2f})\n"
            f"  gap bases {gaps:,} ({100 * gaps / max(c.m_bases, 1):.3f}% of aligned), "
            f"{c.split:,} paths split at a run > max_gap, "
            f"{c.unmodelled:,} I<->D adjacencies\n"
            f"  QV calibration: observed mismatch rate {mm_rate:.2e} against "
            f"mean eps {eps:.2e} from the quality strings\n"
            f"  dropped: {c.skip_mapq:,} mapq, {c.skip_split_read:,} split, "
            f"{c.skip_clip:,} clipped, {c.skip_identity:,} low identity "
            f"({100 * self.kept_frac:.1f}% of primaries kept)"
        )

    @property
    def kept_frac(self) -> float:
        c = self.counts
        seen = c.reads + c.skip_mapq + c.skip_split_read + c.skip_clip + c.skip_identity
        return c.reads / seen if seen else float("nan")


def _ratio(num: int, denom: int) -> tuple[float, float]:
    if denom == 0:
        return 0.0, float("nan")
    p = num / denom
    return p, math.sqrt(max(p * (1 - p), 0.0) / denom)


def estimate(counts: Counts, tau: float = 0.0) -> Estimate:
    """Maximum-likelihood gap parameters for the counted paths.

    `tau` is carried through untouched: read length is conditioned on, so
    termination is not estimated from the data (params.py).
    """
    m_total = counts.m_m + counts.m_i + counts.m_d
    di, se_di = _ratio(counts.m_i, m_total)
    dd, se_dd = _ratio(counts.m_d, m_total)
    gi, se_gi = _ratio(counts.i_i, counts.i_i + counts.i_m)
    gd, se_gd = _ratio(counts.d_d, counts.d_d + counts.d_m)
    return Estimate(PhmmParams(delta_i=di, delta_d=dd, gamma_i=gi, gamma_d=gd, tau=tau),
                    se_di, se_dd, se_gi, se_gd, counts)


def counts_from_paths(paths, max_gap: int = DEFAULT_MAX_GAP) -> Counts:
    """Counts from state strings such as `Sample.path` -- the round-trip route."""
    counts = Counts()
    for path in paths:
        runs, prev = [], None
        for ch in path:
            if ch == prev:
                runs[-1][1] += 1
            else:
                runs.append([ch, 1])
                prev = ch
        op = {"M": _MATCH, "I": _INS, "D": _DEL}
        count_path([(op[s], n) for s, n in runs], counts, max_gap)
        counts.reads += 1
    return counts


# ---- BAM scanning ----------------------------------------------------------

def _windows(bam, n_windows: int, size: int, seed: int):
    """Distinct, non-overlapping windows spread over the genome.

    The genome is tiled into `size`-wide tiles and tiles are drawn *without*
    replacement, for two reasons. Overlapping windows would count a read in the
    intersection twice, and -- worse -- they would break the jackknife, which
    treats each window as an independent block. Tiling also makes the sample
    self-describing: `n_windows` tiles of `size` bp is a stated fraction of the
    assembly, not a number of draws that may or may not land on new sequence.

    Contigs shorter than one tile are skipped, so very short contigs are not
    represented. They are a small share of an assembly's bases and a peculiar
    share of its sequence, which is the trade being made here.
    """
    tiles = [(c, s) for c, L in zip(bam.references, bam.lengths) if L >= size
             for s in range(0, L - size + 1, size)]
    rng = random.Random(seed)
    if n_windows >= len(tiles):
        rng.shuffle(tiles)
        return tiles
    return rng.sample(tiles, n_windows)


def scan_bam(bam_path: str, *, n_windows: int = 500, window: int = 100_000,
             max_reads: int = 20_000, min_mapq: int = 0, max_clip_frac: float = 0.05,
             max_de: float = 0.01, max_gap: int = DEFAULT_MAX_GAP,
             max_per_window: int | None = None, seed: int = 0, progress=None,
             per_window: bool = False):
    """Count gap transitions over primary alignments sampled across the genome.

    The filters keep alignments whose path can be believed: the whole read
    aligned (`max_clip_frac`), no split alignment (an `SA` tag), and high
    identity (`de`, minimap2's gap-compressed divergence).  A read failing those
    is not evidence about sequencing error -- it is evidence about the assembly,
    which is the thing being judged elsewhere.

    `min_mapq` defaults to 0 on purpose.  The target here is a *diploid*
    assembly, where a read from a homozygous region maps equally well to both
    haplotypes and is scored MAPQ 0 however perfect the alignment.  MAPQ
    therefore measures haplotype distinguishability, not correctness, and
    requiring 60 keeps under 10% of primaries -- exactly the divergent ones.
    Measured both ways the gap estimates agree to within ~3 SE, so the filter
    costs precision without buying safety; the mismatch rate, which does differ,
    is what `calibrate_bam` measures properly.

    `max_per_window` caps reads taken from each window, so a global `max_reads`
    spends itself over many windows instead of exhausting the first few. That
    matters because the uncertainty here is *between* windows, not within: reads
    in one window share their sequence context, so they are nowhere near
    independent draws. `per_window=True` returns one `Counts` per window, which
    is what `jackknife` needs to turn that into an honest standard error.
    """
    import pysam

    counts = Counts()
    per: list[Counts] = []
    with pysam.AlignmentFile(bam_path) as bam:
        for contig, start in _windows(bam, n_windows, window, seed):
            here = Counts()
            for a in bam.fetch(contig, start, start + window):
                if counts.reads + here.reads >= max_reads:
                    break
                if max_per_window is not None and here.reads >= max_per_window:
                    break
                if a.is_secondary or a.is_supplementary or a.is_unmapped:
                    continue
                if a.reference_start < start:           # counted in an earlier window
                    continue
                if a.mapping_quality < min_mapq:
                    here.skip_mapq += 1
                    continue
                if a.has_tag("SA"):
                    here.skip_split_read += 1
                    continue
                q_len = a.infer_read_length() or 0
                clipped = q_len - (a.query_alignment_length or 0)
                if not q_len or clipped / q_len > max_clip_frac:
                    here.skip_clip += 1
                    continue
                if a.has_tag("de") and a.get_tag("de") > max_de:
                    here.skip_identity += 1
                    continue

                count_path(a.cigartuples, here, max_gap)
                if a.has_tag("NM"):
                    # NM counts every indel base, including the long runs
                    # count_path drops, so subtract the CIGAR's own gap total
                    # rather than what was counted -- otherwise one 138 bp
                    # deletion reads as 138 mismatches.
                    gaps = sum(n for op, n in a.cigartuples if op in (_INS, _DEL))
                    here.mismatches += max(0, a.get_tag("NM") - gaps)
                if a.query_qualities is not None:
                    eps = eps_from_phred(np.asarray(a.query_qualities, dtype=np.float64))
                    here.eps_sum += float(eps.sum())
                    here.eps_bases += len(eps)
                here.reads += 1
            counts += here
            if here.reads:
                per.append(here)
            if progress is not None and len(per) % 200 == 0:
                progress(counts.reads, max_reads)
            if counts.reads >= max_reads:
                break
    return (counts, per) if per_window else counts


def jackknife(per_window: list) -> dict:
    """Delete-one-window jackknife standard errors for the four parameters.

    The binomial standard error treats every aligned base as an independent
    draw, which they are not: gap rates vary 5x with sequence context, and reads
    in one window share it. Deleting one window at a time and re-estimating
    measures the variation that actually limits the fit, and it is typically an
    order of magnitude larger than the binomial number.
    """
    n = len(per_window)
    if n < 2:
        return {}
    total = Counts()
    for c in per_window:
        total += c
    full = estimate(total).params
    names = ("delta_i", "delta_d", "gamma_i", "gamma_d")
    partials = {k: [] for k in names}
    for skip in per_window:
        left = Counts()
        left += total
        for k, v in vars(skip).items():           # subtract this window
            if isinstance(v, list):
                mine = getattr(left, k)
                for i, x in enumerate(v):
                    mine[i] -= x
            else:
                setattr(left, k, getattr(left, k) - v)
        p = estimate(left).params
        for k in names:
            partials[k].append(getattr(p, k))
    out = {}
    for k in names:
        vals = np.asarray(partials[k])
        out[k] = float(np.sqrt((n - 1) / n * np.sum((vals - vals.mean()) ** 2)))
        out[k + "_estimate"] = float(getattr(full, k))
    out["windows"] = n
    return out


# ---- emission calibration --------------------------------------------------

HOMOPOLYMER_MIN = 4         # a run this long or longer counts as homopolymer context


@dataclass
class Calibration:
    """Observed substitution rate per reported quality, and gap context.

    The pair HMM turns a reported quality straight into an error probability
    (`eps_from_phred`), so a mis-calibrated quality string is a mis-specified
    emission model: every substitution costs `log(3(1-eps)/eps)`, and that
    number is wrong by whatever the calibration is off by.
    """

    bases: np.ndarray = field(default_factory=lambda: np.zeros(MAX_PHRED_BINS, dtype=np.int64))
    mismatches: np.ndarray = field(
        default_factory=lambda: np.zeros(MAX_PHRED_BINS, dtype=np.int64))
    hp_gaps: int = 0            # gap events in a homopolymer of >= HOMOPOLYMER_MIN
    other_gaps: int = 0
    hp_bases: int = 0           # aligned reference bases inside such a homopolymer
    other_bases: int = 0
    reads: int = 0

    def observed_eps(self) -> np.ndarray:
        """Measured substitution rate per reported Phred bin, NaN where unseen."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.bases > 0, self.mismatches / np.maximum(self.bases, 1), np.nan)

    def summary(self) -> str:
        obs = self.observed_eps()
        rows = []
        for q in np.flatnonzero(self.bases):
            if self.bases[q] < 1000:
                continue
            reported = 10.0 ** (-q / 10.0)
            got = obs[q]
            q_obs = -10 * math.log10(got) if got > 0 else float("inf")
            rows.append(f"    Q{q:<3d} {self.bases[q]:>14,}  reported {reported:.2e}  "
                        f"observed {got:.2e}  (Q{q_obs:.1f})")
        hp_rate = self.hp_gaps / max(self.hp_bases, 1)
        other_rate = self.other_gaps / max(self.other_bases, 1)
        return ("quality calibration, {:,} reads\n".format(self.reads)
                + "\n".join(rows)
                + f"\n  gaps in homopolymers >= {HOMOPOLYMER_MIN}: {self.hp_gaps:,} over "
                  f"{self.hp_bases:,} bases ({hp_rate:.2e} per base)\n"
                  f"  gaps elsewhere: {self.other_gaps:,} over {self.other_bases:,} bases "
                  f"({other_rate:.2e} per base) -- "
                  f"{hp_rate / other_rate if other_rate else float('nan'):.1f}x enriched")


def _homopolymer_len(ref: str, i: int) -> int:
    """Length of the homopolymer run of `ref` containing position i."""
    if not 0 <= i < len(ref):
        return 0
    base, lo, hi = ref[i], i, i
    while lo > 0 and ref[lo - 1] == base:
        lo -= 1
    while hi + 1 < len(ref) and ref[hi + 1] == base:
        hi += 1
    return hi - lo + 1


def calibrate_bam(bam_path: str, fasta_path: str, *, n_windows: int = 100,
                  window: int = 100_000, max_reads: int = 2_000, min_mapq: int = 0,
                  max_de: float = 0.01, max_gap: int = DEFAULT_MAX_GAP,
                  seed: int = 0, progress=None) -> Calibration:
    """Compare read to reference base by base, tallying by reported quality.

    Walks the CIGAR directly rather than using `get_aligned_pairs(with_seq=True)`,
    which needs an MD tag these BAMs do not carry.
    """
    import pysam

    cal = Calibration()
    with pysam.AlignmentFile(bam_path) as bam, pysam.FastaFile(fasta_path) as fa:
        for contig, start in _windows(bam, n_windows, window, seed):
            for a in bam.fetch(contig, start, start + window):
                if cal.reads >= max_reads:
                    return cal
                if (a.is_secondary or a.is_supplementary or a.is_unmapped
                        or a.reference_start < start or a.mapping_quality < min_mapq
                        or a.has_tag("SA") or a.query_qualities is None):
                    continue
                if a.has_tag("de") and a.get_tag("de") > max_de:
                    continue
                ref = fa.fetch(a.reference_name, a.reference_start, a.reference_end).upper()
                read = a.query_sequence.upper()
                qual = np.asarray(a.query_qualities, dtype=np.int64)
                ri, qi = 0, 0                    # ref offset, read offset
                for op, n in a.cigartuples:
                    if op in (_MATCH, _EQ, _DIFF):
                        r = np.frombuffer(ref[ri:ri + n].encode(), dtype=np.uint8)
                        q = np.frombuffer(read[qi:qi + n].encode(), dtype=np.uint8)
                        qv = qual[qi:qi + n]
                        np.add.at(cal.bases, qv, 1)
                        np.add.at(cal.mismatches, qv[r != q], 1)
                        ri += n
                        qi += n
                    elif op == _INS:
                        if n <= max_gap:
                            hp = _homopolymer_len(ref, ri - 1) >= HOMOPOLYMER_MIN
                            cal.hp_gaps += hp
                            cal.other_gaps += not hp
                        qi += n
                    elif op == _DEL:
                        if n <= max_gap:
                            hp = _homopolymer_len(ref, ri) >= HOMOPOLYMER_MIN
                            cal.hp_gaps += hp
                            cal.other_gaps += not hp
                        ri += n
                    elif op in (_SOFT,):
                        qi += n
                    elif op == _REF_SKIP:
                        ri += n
                # Homopolymer share of the aligned reference, the denominator for
                # the gap rates above.
                arr = np.frombuffer(ref.encode(), dtype=np.uint8)
                if len(arr) > 1:
                    same = arr[1:] == arr[:-1]
                    run, hp = 1, 0
                    for s in same:
                        if s:
                            run += 1
                        else:
                            hp += run if run >= HOMOPOLYMER_MIN else 0
                            run = 1
                    hp += run if run >= HOMOPOLYMER_MIN else 0
                    cal.hp_bases += hp
                    cal.other_bases += len(arr) - hp
                cal.reads += 1
                if progress is not None and cal.reads % 500 == 0:
                    progress(cal.reads, max_reads)
    return cal


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description="Estimate pair-HMM gap parameters from a BAM")
    ap.add_argument("bam")
    ap.add_argument("--out", help="write the estimate as JSON")
    ap.add_argument("--windows", type=int, default=500)
    ap.add_argument("--window-size", type=int, default=100_000)
    ap.add_argument("--max-reads", type=int, default=20_000)
    ap.add_argument("--max-per-window", type=int, default=None,
                    help="cap reads taken from each window, so the budget spreads")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--max-de", type=float, default=0.01)
    ap.add_argument("--max-gap", type=int, default=DEFAULT_MAX_GAP)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--calibrate", metavar="FASTA",
                    help="also measure substitution rate per reported quality "
                         "against this reference")
    ap.add_argument("--calibrate-reads", type=int, default=2000)
    args = ap.parse_args(argv)

    counts, per = scan_bam(args.bam, n_windows=args.windows, window=args.window_size,
                           max_reads=args.max_reads, min_mapq=args.min_mapq,
                           max_de=args.max_de, max_gap=args.max_gap, seed=args.seed,
                           max_per_window=args.max_per_window, per_window=True,
                           progress=lambda i, n: print(f"  {i:,}/{n:,} reads", flush=True))
    est = estimate(counts)
    print(est.summary())
    jk = jackknife(per)
    if jk:
        print(f"  jackknife over {jk['windows']:,} windows: "
              f"delta_i +/- {jk['delta_i']:.1e}, delta_d +/- {jk['delta_d']:.1e}, "
              f"gamma_i +/- {jk['gamma_i']:.4f}, gamma_d +/- {jk['gamma_d']:.4f}")
    cal = None
    if args.calibrate:
        cal = calibrate_bam(args.bam, args.calibrate, n_windows=args.windows,
                            window=args.window_size, max_reads=args.calibrate_reads,
                            min_mapq=args.min_mapq, max_de=args.max_de,
                            max_gap=args.max_gap, seed=args.seed)
        print(cal.summary())
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"params": asdict(est.params),
                       "se": {"delta_i": est.se_delta_i, "delta_d": est.se_delta_d,
                              "gamma_i": est.se_gamma_i, "gamma_d": est.se_gamma_d},
                       "counts": vars(counts),
                       "jackknife": jk,
                       "calibration": (None if cal is None else
                                       {"bases": cal.bases.tolist(),
                                        "mismatches": cal.mismatches.tolist(),
                                        "hp_gaps": cal.hp_gaps, "other_gaps": cal.other_gaps,
                                        "hp_bases": cal.hp_bases,
                                        "other_bases": cal.other_bases,
                                        "reads": cal.reads}),
                       "filters": {"bam": args.bam, "min_mapq": args.min_mapq,
                                   "max_de": args.max_de, "max_gap": args.max_gap,
                                   "windows": args.windows, "seed": args.seed}},
                      fh, indent=1)
        print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
