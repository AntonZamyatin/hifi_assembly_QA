"""Tests for the row-recycled forward and the vectorised deletion row (F1).

Both changes are pure optimisation: they must return exactly what the
full-table reference returns. So the tests here are equivalence tests against
`_dp`, plus the scale check that motivated the work -- the reference cannot
allocate a production-size read at all, and this path can.

Run: venv/bin/python tests/test_phmm_rows.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.phmm import PhmmParams, quality_table
from qamodel.phmm.dp import (NEG_INF, _deletion_row, _dp_two_rows, _prepare,
                             _restrict_to_end, _score_from_terminal, forward)


def forward_rows(read, qual, template, params, qt, *, origin=None,
                 full_template=False):
    """The two-row engine, wrapped for comparison against `forward`.

    `forward_window` deliberately omits `full_template`, so this reaches the
    engine directly -- the point here is to check the engine in both modes, not
    to exercise the public wrapper.
    """
    r, q, t, p, tab = _prepare(read, qual, template, params, qt)
    terminal, _ = _dp_two_rows(r, q, t, p, tab, origin, None)
    if full_template:
        terminal = _restrict_to_end(terminal, t.shape[0])
    return type("R", (), {"logp": _score_from_terminal(terminal, False),
                          "terminal": terminal})

QT = quality_table()
_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _scalar_deletion_row(cur_m, m_d, d_d):
    """The recurrence written out, one cell at a time. Ground truth."""
    out = np.full(cur_m.shape[0], NEG_INF)
    prev = NEG_INF
    for v in range(1, cur_m.shape[0]):
        prev = np.logaddexp(cur_m[v - 1] + m_d, prev + d_d)
        out[v] = prev
    return out


def test_deletion_row_matches_the_scalar_recurrence():
    """The unrolled form must reproduce the left-to-right scan exactly."""
    rng = np.random.default_rng(4)
    worst = 0.0
    for n in (2, 5, 50, 400, 2000):
        for gd, dd in ((0.3, 1e-4), (0.9, 1e-2), (0.01, 0.5)):
            p = PhmmParams(delta_d=dd, gamma_d=gd)
            tr = p.log_trans()
            cur = np.log(rng.random(n))
            cur[rng.random(n) < 0.2] = NEG_INF          # dead cells are normal
            want = _scalar_deletion_row(cur, tr.m_d, tr.d_d)
            got = _deletion_row(cur, tr.m_d, tr.d_d)
            fin = np.isfinite(want)
            d = np.abs(want[fin] - got[fin]).max() if fin.any() else 0.0
            worst = max(worst, d)
            check(d < 1e-9, f"n={n} gd={gd}: deletion row differs by {d:.2e}")
            check(np.array_equal(np.isfinite(want), np.isfinite(got)),
                  f"n={n} gd={gd}: -inf pattern differs")
    print(f"    worst |delta| over all sizes and gap settings = {worst:.2e}")


def test_deletion_row_degenerate_parameters():
    """gamma_D = 0 means gaps cannot extend; delta_D = 0 means none can open."""
    cur = np.array([-1.0, -2.0, -3.0, -4.0])
    p = PhmmParams(delta_d=1e-3, gamma_d=0.0)
    tr = p.log_trans()
    got = _deletion_row(cur, tr.m_d, tr.d_d)
    check(np.allclose(got[1:], cur[:-1] + tr.m_d),
          f"gamma_D = 0 must give length-1 deletions only, got {got}")
    p0 = PhmmParams(delta_d=0.0)
    got0 = _deletion_row(cur, p0.log_trans().m_d, p0.log_trans().d_d)
    check(np.all(np.isneginf(got0)), f"delta_D = 0 must forbid every D, got {got0}")


def test_matches_reference_across_configurations():
    """The anchor for F1: same number as the full-table path, everywhere."""
    rng = random.Random(99)
    settings = [
        PhmmParams(),
        PhmmParams(delta_i=0.02, delta_d=0.03, gamma_i=0.25, gamma_d=0.35),
        PhmmParams(delta_i=0.0, delta_d=0.0),
        PhmmParams(delta_d=1e-3, gamma_d=0.0),
        PhmmParams(tau=1e-3),
    ]
    cases = 0
    worst = 0.0
    for params in settings:
        for _ in range(6):
            L = rng.randint(1, 25)
            T = rng.randint(1, 40)
            read, tmpl = _rand_seq(L, rng), _rand_seq(T, rng)
            qual = [rng.choice((15, 30, 40)) for _ in range(L)]
            origins = [None, 0, min(T - 1, rng.randint(0, max(0, T - 1)))]
            if T >= 3:
                origins.append((1, T - 1))
            for origin in origins:
                for full in (False, True):
                    a = forward(read, qual, tmpl, params, QT,
                                origin=origin, full_template=full)
                    b = forward_rows(read, qual, tmpl, params, QT,
                                     origin=origin, full_template=full)
                    cases += 1
                    if not math.isfinite(a.logp) or not math.isfinite(b.logp):
                        check(a.logp == b.logp,
                              f"impossible case disagrees: {a.logp} vs {b.logp}")
                        continue
                    d = abs(a.logp - b.logp)
                    worst = max(worst, d)
                    check(d < 1e-9,
                          f"logp differs by {d:.2e} (L={L} T={T} origin={origin} "
                          f"full={full})")
                    td = np.abs(a.terminal[np.isfinite(a.terminal)]
                                - b.terminal[np.isfinite(a.terminal)])
                    check(td.max() < 1e-9 if td.size else True,
                          "terminal distribution differs")
    print(f"    {cases} configurations, worst |delta| = {worst:.2e}")


def test_matches_reference_at_moderate_scale():
    """Small cases can hide accumulation error; this one would not."""
    rng = random.Random(2)
    tmpl = _rand_seq(3000, rng)
    read = tmpl[500:1200]
    read = read[:300] + "A" + read[300:450] + read[470:]      # an indel or two
    a = forward(read, 38, tmpl, PhmmParams(), QT).logp
    b = forward_rows(read, 38, tmpl, PhmmParams(), QT).logp
    check(abs(a - b) < 1e-9, f"at L={len(read)} T=3000 the two differ by {abs(a-b):.2e}")
    print(f"    L={len(read)} T=3000: reference {a:.6f}, rows {b:.6f}")


def test_runs_at_a_size_the_reference_cannot():
    """The point of the work: 5.1 GB of matrices versus two rows."""
    rng = random.Random(17)
    tmpl = _rand_seq(12_000, rng)
    read = tmpl[100:11_900]
    t0 = time.perf_counter()
    got = forward_rows(read, 40, tmpl, PhmmParams(), QT)
    dt = time.perf_counter() - t0
    full_gb = 3 * 8 * (len(read) + 1) * (len(tmpl) + 1) / 2**30
    check(math.isfinite(got.logp), f"score must be finite, got {got.logp}")
    check(got.logp > -50.0, f"a near-exact read should score near zero, got {got.logp}")
    print(f"    L={len(read)} T={len(tmpl)}: logp={got.logp:.4f} in {dt:.1f}s "
          f"(full matrices would need {full_gb:.1f} GB)")


def main():
    tests = [
        test_deletion_row_matches_the_scalar_recurrence,
        test_deletion_row_degenerate_parameters,
        test_matches_reference_across_configurations,
        test_matches_reference_at_moderate_scale,
        test_runs_at_a_size_the_reference_cannot,
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
