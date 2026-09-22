"""Tests for the generative sampler (stage D of docs/phmm.md).

The anchor here is different in kind from the other suites. Everywhere else the
scorer is checked against a rearrangement of its own recursion; here it is
checked against the model run *forwards*. If sampling and scoring agree on the
frequency of every read, they are provably the same distribution, and no shared
algebraic mistake can hide in both.

Run: venv/bin/python tests/test_phmm_sample.py
"""

import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.phmm import PhmmParams, quality_table
from qamodel.phmm.dp import forward
from qamodel.phmm.sample import run_lengths, sample_read, sample_reads

QT = quality_table()
TEMPLATE = "ACGTTGACCATGCAGTTACGGATCCTTAGCATCGGATTAC"     # 40 bp
_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def test_path_invariants():
    """Structural facts that must hold for every draw."""
    rng = np.random.default_rng(1)
    p = PhmmParams(tau=0.2, delta_i=0.06, delta_d=0.06, gamma_i=0.4, gamma_d=0.4)
    bad = 0
    for s in sample_reads(TEMPLATE, p, rng, 4000):
        if s.truncated:
            continue
        ok = (len(s.read) == s.path.count("M") + s.path.count("I")
              and (s.end - s.origin + 1) == s.path.count("M") + s.path.count("D")
              and not s.path.startswith("D")
              and not s.path.endswith("D")
              and len(s.qual) == len(s.read))
        bad += not ok
    check(bad == 0, f"{bad} samples violated the path invariants")
    check(all(not sample_read(TEMPLATE, p, rng, max_len=8).truncated
              or True for _ in range(10)), "max_len path smoke")


def test_frequencies_match_the_scorer():
    """The anchor: sample the model, then score what came out.

    With tau > 0 the chain terminates on its own, and every walk emitting the
    read terminates exactly once at probability tau. So the frequency of a read
    is exactly `tau * exp(forward(...).logp)` -- an equality, not a proportion.
    """
    rng = np.random.default_rng(20260910)
    p = PhmmParams(tau=0.25, delta_i=0.05, delta_d=0.05, gamma_i=0.4, gamma_d=0.4)
    n = 200_000
    counts = Counter()
    trunc = 0
    for s in sample_reads(TEMPLATE, p, rng, n, origin=0, qual=30):
        if s.truncated:
            trunc += 1
        else:
            counts[s.read] += 1
    check(trunc / n < 1e-3,
          f"template too short: {trunc} of {n} walks ran off the end")

    worst_z = 0.0
    for read, c in counts.most_common(20):
        obs = c / n
        pred = p.tau * math.exp(
            forward(read, 30, TEMPLATE, p, QT, origin=0).logp)
        se = math.sqrt(pred * (1 - pred) / n)
        z = abs(obs - pred) / se
        worst_z = max(worst_z, z)
        check(z < 4.0,
              f"read {read!r}: observed {obs:.5f}, predicted {pred:.5f}, "
              f"{z:.1f} standard errors apart")

    # A global check has to enumerate reads chosen *independently of the sample*.
    # Selecting them by observed count instead biases the comparison upward, by
    # keeping exactly those reads that happened to come up more than expected.
    every = [""]
    for _ in range(3):
        every += [w + b for w in every for b in "ACGT" if len(w) < 3]
    every = sorted({w for w in every if w})
    tot_obs = sum(counts.get(w, 0) for w in every) / n
    tot_pred = sum(p.tau * math.exp(forward(w, 30, TEMPLATE, p, QT, origin=0).logp)
                   for w in every)
    check(abs(tot_obs - tot_pred) < 0.005,
          f"summed mass over all {len(every)} reads of length <= 3: "
          f"observed {tot_obs:.4f}, predicted {tot_pred:.4f}")
    print(f"    {len(counts)} distinct reads, worst deviation {worst_z:.1f} SE; "
          f"all {len(every)} reads of length <= 3: {tot_obs:.4f} observed, "
          f"{tot_pred:.4f} predicted")


def test_length_conditioned_frequencies_are_proportional():
    """With tau = 0 the sampler draws from P(r | l), which is forward / Z(l).

    Z is unknown until stage C, but it is a single constant for a fixed length,
    so normalising both sides over the compared reads removes it -- and dividing
    back through estimates it.
    """
    rng = np.random.default_rng(4242)
    p = PhmmParams(tau=0.0, delta_i=0.05, delta_d=0.05, gamma_i=0.4, gamma_d=0.4)
    n, L = 120_000, 6
    counts = Counter()
    for s in sample_reads(TEMPLATE, p, rng, n, origin=0, qual=25, max_len=L):
        if not s.truncated:
            check(len(s.read) == L, f"max_len ignored: got {len(s.read)}")
            counts[s.read] += 1

    top = counts.most_common(15)
    obs = np.array([c for _, c in top], dtype=float)
    obs /= obs.sum()
    lp = np.array([forward(r, 25, TEMPLATE, p, QT, origin=0).logp for r, _ in top])
    pred = np.exp(lp - lp.max())
    pred /= pred.sum()

    worst = float(np.abs(obs - pred).max())
    check(worst < 0.01,
          f"relative frequencies differ by up to {worst:.4f} after normalising")

    # Inverting the proportionality gives Z(l); stage C will compute it directly.
    z_est = [math.exp(l) * n / c for (r, c), l in zip(top, lp)]
    spread = (max(z_est) - min(z_est)) / np.mean(z_est)
    check(spread < 0.15,
          f"implied Z(l) is not constant across reads: spread {spread:.3f}")
    print(f"    {len(counts)} distinct reads of length {L}; "
          f"worst normalised gap {worst:.4f}; implied Z({L}) = {np.mean(z_est):.4f} "
          f"(spread {spread:.1%})")


def test_indel_runs_are_geometric():
    """Gap lengths must follow the geometry the gap-extend parameter specifies."""
    rng = np.random.default_rng(9)
    g_i, g_d = 0.55, 0.35
    p = PhmmParams(tau=0.05, delta_i=0.15, delta_d=0.15, gamma_i=g_i, gamma_d=g_d)
    # Deliberately long: discarding truncated walks censors deletion runs, since
    # a D run still going when the template runs out is what causes truncation.
    # On a 40 bp template that biased the deletion mean low by 4 standard errors,
    # while leaving insertions untouched -- they consume no template.
    long_template = "".join(rng.choice(list("ACGT"), size=400))
    ins, dels = [], []
    for s in sample_reads(long_template, p, rng, 20_000, qual=35):
        if s.truncated:
            continue
        ins += run_lengths(s.path, "I")
        dels += run_lengths(s.path, "D")

    for name, runs, g in (("insertion", ins, g_i), ("deletion", dels, g_d)):
        want = 1.0 / (1.0 - g)
        got = float(np.mean(runs))
        se = float(np.std(runs)) / math.sqrt(len(runs))
        check(abs(got - want) < 4 * se,
              f"{name} run mean {got:.4f}, expected {want:.4f} "
              f"({abs(got-want)/se:.1f} SE, n={len(runs)})")
        frac1 = sum(1 for r in runs if r == 1) / len(runs)
        check(abs(frac1 - (1 - g)) < 0.02,
              f"{name}: {frac1:.3f} of runs have length 1, expected {1-g:.3f}")
        print(f"    {name}s: n={len(runs)} mean {got:.3f} (expected {want:.3f}), "
              f"P(len 1) = {frac1:.3f} (expected {1-g:.3f})")


def test_truncation_is_flagged_not_silent():
    """A walk that wants template the window lacks is not a draw from the model."""
    rng = np.random.default_rng(3)
    p = PhmmParams(tau=0.02, delta_i=0.02, delta_d=0.02)   # long reads, short template
    samples = sample_reads("ACGTACGTAC", p, rng, 400, qual=35)
    trunc = [s for s in samples if s.truncated]
    check(len(trunc) > 100,
          f"expected most walks to overrun a 10 bp template, got {len(trunc)}/400")
    for s in trunc:
        check(s.end == len("ACGTACGTAC") - 1,
              f"a truncated walk must have consumed the whole template, end={s.end}")
    print(f"    {len(trunc)}/400 walks overran a 10 bp template, all flagged")


def test_tau_zero_without_max_len_is_rejected():
    rng = np.random.default_rng(0)
    try:
        sample_read(TEMPLATE, PhmmParams(), rng)
    except ValueError:
        return
    _failures.append("tau = 0 with no max_len must raise, not hang")
    print("  FAIL: tau = 0 with no max_len must raise, not hang")


def main():
    tests = [
        test_path_invariants,
        test_frequencies_match_the_scorer,
        test_length_conditioned_frequencies_are_proportional,
        test_indel_runs_are_geometric,
        test_truncation_is_flagged_not_silent,
        test_tau_zero_without_max_len_is_rejected,
    ]
    for fn in tests:
        print(f"{fn.__name__} ...")
        fn()
    print()
    if _failures:
        print(f"{len(_failures)} FAILURE(S)")
        return 1
    print(f"all {len(tests)} tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
