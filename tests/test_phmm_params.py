"""Tests for pair HMM parameters and emissions (Stage A of docs/phmm.md).

The substantive test here is that the exponentiated log transitions form a
proper stochastic matrix.  That checks two things at once: the table of
model.tex §5.3 is a distribution at all, and the linear -> log conversion is
faithful.  Everything downstream is a sum over paths through this matrix, so an
error here would show up only as a plausible-looking wrong number.

Run: venv/bin/python tests/test_phmm_params.py
"""

import math
import os
import sys
from dataclasses import fields

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.phmm import (
    EPS_MAX,
    MAX_PHRED,
    LogTrans,
    PhmmParams,
    eps_from_phred,
    mismatch_cost,
    quality_table,
)

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def test_mismatch_cost_matches_doc():
    """model.tex §5.4: 5.7 nats at Q20, 8.0 at Q30, 10.3 at Q40."""
    for q, want in ((20, 5.7), (30, 8.0), (40, 10.3)):
        got = float(mismatch_cost(q))
        check(abs(got - want) < 0.05, f"s(Q{q}) = {got:.3f}, expected {want}")

    # The placeholder being replaced: exp(-NM) charges 1 nat per mismatch.
    implied = 3.0 / (3.0 + math.e)
    check(abs(float(mismatch_cost(-10 * math.log10(implied))) - 1.0) < 1e-9,
          "s = 1 must correspond to eps = 3/(3+e)")


def test_transition_rows_are_distributions():
    p = PhmmParams()
    t = p.log_trans()
    rows = {
        "B": (t.b_m, t.b_i),
        "M": (t.m_m, t.m_i, t.m_d, t.m_e),
        "I": (t.i_i, t.i_m, t.i_e),
        "D": (t.d_d, t.d_m),
    }
    for name, logs in rows.items():
        total = sum(math.exp(x) for x in logs)
        check(abs(total - 1.0) < 1e-12,
              f"row {name} sums to {total!r}, must be 1")


def test_forbidden_transitions_absent():
    """B->D and D->E are excluded by §5.3; the type must not offer them."""
    names = {f.name for f in fields(LogTrans)}
    check("d_e" not in names, "LogTrans must not expose D->E (termination from D)")
    check("b_d" not in names, "LogTrans must not expose B->D (duplicates origin s+1)")
    check(len(names) == 11, f"expected 11 edges, got {len(names)}: {sorted(names)}")


def test_quality_table_agrees_with_mismatch_cost():
    """Two routes to s must not drift apart."""
    qt = quality_table()
    q = np.arange(MAX_PHRED + 1)
    check(np.allclose(qt.cost(), mismatch_cost(q), atol=1e-12),
          "QualityTable.cost() must equal mismatch_cost()")
    check(abs(qt.log_insert - math.log(0.25)) < 1e-12,
          "default insertion emission must be log(1/4)")


def test_emissions_are_finite_and_normalised():
    """Match plus three mismatches must sum to one at every quality."""
    qt = quality_table()
    total = np.exp(qt.log_match) + 3.0 * np.exp(qt.log_mismatch)
    check(np.allclose(total, 1.0, atol=1e-12),
          f"emissions must sum to 1; worst deviation {np.abs(total - 1).max():.2e}")
    check(np.all(np.isfinite(qt.log_match)) and np.all(np.isfinite(qt.log_mismatch)),
          "no emission may be -inf or nan, including at Q0")


def test_quality_clamping():
    """Q0 means 'unknown' in HiFi BAMs; unclamped it gives eps = 1 and log(0)."""
    check(float(eps_from_phred(0)) == EPS_MAX,
          f"Q0 must clamp to EPS_MAX, got {float(eps_from_phred(0))}")
    check(float(mismatch_cost(0)) == 0.0,
          "at the clamp a mismatch must cost exactly zero, not go negative")
    check(np.all(np.diff(mismatch_cost(np.arange(MAX_PHRED + 1))) >= 0),
          "s must be non-decreasing in reported quality")
    check(np.isfinite(mismatch_cost(MAX_PHRED)),
          "s must stay finite at the top of the Phred range")


def test_validation_rejects_impossible_parameters():
    bad = [
        (dict(delta_i=0.6, delta_d=0.6), "M->M non-positive"),
        (dict(gamma_i=0.999, tau=0.01), "I->M non-positive"),
        (dict(gamma_d=1.0), "gamma_d = 1 leaves D with no exit"),
        (dict(tau=-0.1), "negative probability"),
    ]
    for kwargs, why in bad:
        try:
            PhmmParams(**kwargs)
        except ValueError:
            continue
        _failures.append(f"PhmmParams({kwargs}) must raise: {why}")
        print(f"  FAIL: PhmmParams({kwargs}) must raise: {why}")

    PhmmParams(delta_i=0.0, delta_d=0.0)      # gapless is degenerate but legal


def main():
    tests = [
        test_mismatch_cost_matches_doc,
        test_transition_rows_are_distributions,
        test_forbidden_transitions_absent,
        test_quality_table_agrees_with_mismatch_cost,
        test_emissions_are_finite_and_normalised,
        test_quality_clamping,
        test_validation_rejects_impossible_parameters,
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
