"""Tests for gap-parameter estimation (stage G of docs/phmm.md).

The round trip is the real test: reads sampled by stage D carry the path that
produced them, so estimating from those paths must recover the parameters that
generated them. That is an identifiability check, and it is the one place where
the path being counted really is the model's own -- a CIGAR's is the aligner's.

The rest pin the counting rules that keep a real BAM from poisoning the fit:
runs too long to be sequencing indels split the path, clips break it, and the
I<->D adjacencies the model has no edge for are counted, not folded away.

Run: venv/bin/python tests/test_estimate.py
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.estimate import (
    Counts,
    count_path,
    counts_from_paths,
    estimate,
)
from qamodel.phmm import PhmmParams, sample_reads

_failures = []

_OPS = {"M": 0, "I": 1, "D": 2, "S": 4, "H": 5, "N": 3}


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def cigar(text: str):
    """"5M2I3M" -> pysam-style (op, length) tuples."""
    out, num = [], ""
    for ch in text:
        if ch.isdigit():
            num += ch
        else:
            out.append((_OPS[ch], int(num)))
            num = ""
    return out


def test_transition_counts():
    c = Counts()
    count_path(cigar("5M2I3M1D4M"), c)
    check((c.m_bases, c.i_bases, c.d_bases) == (12, 2, 1),
          f"bases per state: {(c.m_bases, c.i_bases, c.d_bases)}")
    check(c.m_m == 4 + 2 + 3 and c.i_i == 1 and c.d_d == 0,
          f"within-run transitions: m_m {c.m_m}, i_i {c.i_i}, d_d {c.d_d}")
    check((c.m_i, c.i_m, c.m_d, c.d_m) == (1, 1, 1, 1),
          f"between-run transitions: {(c.m_i, c.i_m, c.m_d, c.d_m)}")
    check(c.i_hist[2] == 1 and c.d_hist[1] == 1, "run lengths tabulated")


def test_estimator_is_the_transition_ratio():
    c = Counts()
    count_path(cigar("100M1I100M2D100M"), c)
    e = estimate(c)
    m_total = c.m_m + c.m_i + c.m_d
    check(abs(e.params.delta_i - 1 / m_total) < 1e-12, "delta_i = n(M->I)/n(M->.)")
    check(abs(e.params.delta_d - 1 / m_total) < 1e-12, "delta_d = n(M->D)/n(M->.)")
    check(e.params.gamma_i == 0.0, "a 1 bp insertion never extends")
    check(abs(e.params.gamma_d - 0.5) < 1e-12, "a 2 bp deletion extends once in two")
    check(e.params.tau == 0.0, "tau is conditioned on, never estimated")


def test_long_runs_split_the_path():
    """A 138 bp deletion is a variant, not a sequencing indel."""
    c = Counts()
    count_path(cigar("100M138D100M"), c, max_gap=20)
    check(c.d_bases == 0 and c.split == 1, f"the long run is dropped: {c.d_bases}, {c.split}")
    check(c.m_d == 0 and c.d_m == 0, "and so are the transitions touching it")
    check(c.m_m == 99 + 99, f"both sides still counted as paths: m_m {c.m_m}")
    check(estimate(c).params.gamma_d == 0.0, "gamma_d is not inflated towards 1")

    kept = Counts()
    count_path(cigar("100M3D100M"), kept, max_gap=20)
    check(kept.d_bases == 3 and kept.split == 0, "a short run is kept")


def test_clips_break_the_path():
    c = Counts()
    count_path(cigar("10S5M3S"), c)
    check(c.m_bases == 5 and c.m_m == 4, "clipped bases are not aligned bases")
    check(c.m_i == c.i_m == c.m_d == c.d_m == 0,
          "a clip is not a transition: the free start and end already handle it")


def test_unmodelled_adjacency_is_counted_not_folded():
    c = Counts()
    count_path(cigar("10M2I3D10M"), c)
    check(c.i_d == 1 and c.unmodelled == 1, f"I->D counted: {c.i_d}")
    check(c.d_i == 0 and c.i_m == 0, "and not silently turned into I->M")


def test_round_trip_recovers_generating_parameters():
    """Estimating from the sampler's own paths must return what generated them."""
    truth = PhmmParams(delta_i=0.004, delta_d=0.006, gamma_i=0.25, gamma_d=0.40)
    rng = np.random.default_rng(7)
    tmpl = "".join(rng.choice(list("ACGT"), 20_000))
    reads = sample_reads(tmpl, truth, rng, 200, origin=0, max_len=5_000)
    paths = [s.path for s in reads if not s.truncated]
    check(len(paths) > 150, f"enough untruncated samples: {len(paths)}")

    e = estimate(counts_from_paths(paths))
    for name, se in (("delta_i", e.se_delta_i), ("delta_d", e.se_delta_d),
                     ("gamma_i", e.se_gamma_i), ("gamma_d", e.se_gamma_d)):
        got, want = getattr(e.params, name), getattr(truth, name)
        z = abs(got - want) / se
        check(z < 4.0, f"{name}: {got:.5f} vs {want} is {z:.1f} SE out (se {se:.5f})")


def main():
    tests = [test_transition_counts, test_estimator_is_the_transition_ratio,
             test_long_runs_split_the_path, test_clips_break_the_path,
             test_unmodelled_adjacency_is_counted_not_folded,
             test_round_trip_recovers_generating_parameters]
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
