"""Tests for the compiled engine (stage F4 of docs/phmm.md).

The kernel computes the same recursion in linear space with per-row rescaling
rather than in log space, so it is a genuinely different implementation and not
a transcription. That makes agreement with `dp.py` the whole test: every check
here is "does it return what the oracle returns", including in the degenerate
parameter corners where the rescaling could divide by zero.

Run: venv/bin/python tests/test_phmm_kernel.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.phmm import PhmmParams, quality_table
from qamodel.phmm.band import Band
from qamodel.phmm.dp import forward, forward_window
from qamodel.phmm.kernel import HAVE_NUMBA, forward_numba

QT = quality_table()
_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_numba_is_available():
    check(HAVE_NUMBA, "numba missing: the kernel falls back to pure python and "
                      "the timings below are meaningless")


def test_matches_the_oracle():
    """Against the full-table log-space engine, over many small configurations."""
    rng = random.Random(7)
    settings = [
        PhmmParams(),
        PhmmParams(delta_i=0.02, delta_d=0.03, gamma_i=0.25, gamma_d=0.35),
        PhmmParams(delta_i=0.0, delta_d=0.0),          # gapless
        PhmmParams(delta_d=1e-3, gamma_d=0.0),         # gaps cannot extend
        PhmmParams(delta_i=0.0),                       # no insertions
        PhmmParams(tau=1e-3),                          # joint model
    ]
    cases, worst = 0, 0.0
    for params in settings:
        for _ in range(6):
            T = rng.randint(1, 40)
            L = rng.randint(1, 25)
            read, tmpl = _rand_seq(L, rng), _rand_seq(T, rng)
            qual = [rng.choice((12, 30, 40)) for _ in range(L)]
            origins = [None, 0]
            if T >= 3:
                origins += [(1, T - 1), T - 1]
            for origin in origins:
                want = forward(read, qual, tmpl, params, QT, origin=origin).logp
                got = forward_numba(read, qual, tmpl, params, QT, origin=origin).logp
                cases += 1
                if not math.isfinite(want) or not math.isfinite(got):
                    check(want == got, f"impossible case disagrees: {want} vs {got}")
                    continue
                d = abs(want - got)
                worst = max(worst, d)
                check(d < 1e-9,
                      f"logp differs by {d:.2e} (L={L} T={T} origin={origin} "
                      f"params={params})")
    print(f"    {cases} configurations, worst |delta| = {worst:.2e}")


def test_terminal_row_matches():
    """Not just the total: the whole distribution over end positions."""
    rng = random.Random(31)
    tmpl = _rand_seq(300, rng)
    read = tmpl[60:200]
    a = forward_window(read, 37, tmpl, PhmmParams(), QT)
    b = forward_numba(read, 37, tmpl, PhmmParams(), QT)
    fin = np.isfinite(a.terminal) & np.isfinite(b.terminal)
    d = np.abs(a.terminal[fin] - b.terminal[fin]).max()
    check(d < 1e-9, f"terminal rows differ by {d:.2e}")
    check(abs(a.edge_mass(2) - b.edge_mass(2)) < 1e-9,
          "edge_mass must agree")
    print(f"    worst |delta| across {fin.sum()} end positions = {d:.2e}")


def test_banded_matches():
    rng = random.Random(64)
    T = 1200
    tmpl = _rand_seq(T, rng)
    read = tmpl[100:900]
    L = len(read)
    for width in (20, 100, 400):
        band = Band.diagonal(L, T, 100, 101, width=width)
        a = forward_window(read, 40, tmpl, PhmmParams(), QT,
                           origin=(100, 101), band=band)
        b = forward_numba(read, 40, tmpl, PhmmParams(), QT,
                          origin=(100, 101), band=band)
        check(abs(a.logp - b.logp) < 1e-9,
              f"width {width}: logp differs by {abs(a.logp-b.logp):.2e}")
        check(abs(a.band_edge - b.band_edge) < 1e-9,
              f"width {width}: band_edge differs "
              f"({a.band_edge:.3e} vs {b.band_edge:.3e})")
    print("    diagonal bands at widths 20, 100, 400 all agree")


def test_stable_on_a_long_read():
    """Rescaling must hold over thousands of rows, where log space cannot drift."""
    rng = random.Random(5)
    T = 9000
    tmpl = _rand_seq(T, rng)
    read = tmpl[50:8800]
    a = forward_window(read, 40, tmpl, PhmmParams(), QT, origin=(50, 51))
    b = forward_numba(read, 40, tmpl, PhmmParams(), QT, origin=(50, 51))
    check(abs(a.logp - b.logp) < 1e-6,
          f"L={len(read)}: {a.logp:.6f} vs {b.logp:.6f}, delta {abs(a.logp-b.logp):.2e}")
    print(f"    L={len(read)} T={T}: numpy {a.logp:.6f}, numba {b.logp:.6f}")


def test_speedup():
    rng = random.Random(11)
    T = 6000
    tmpl = _rand_seq(T, rng)
    read = tmpl[100:5900]
    L = len(read)
    band = Band.diagonal(L, T, 100, 101, width=100)
    kw = dict(origin=(100, 101), band=band)

    forward_numba(read[:50], 40, tmpl[:200], PhmmParams(), QT)   # warm the JIT

    t0 = time.perf_counter()
    a = forward_window(read, 40, tmpl, PhmmParams(), QT, **kw)
    t1 = time.perf_counter()
    b = forward_numba(read, 40, tmpl, PhmmParams(), QT, **kw)
    t2 = time.perf_counter()

    check(abs(a.logp - b.logp) < 1e-9, f"answers differ by {abs(a.logp-b.logp):.2e}")
    speedup = (t1 - t0) / (t2 - t1)
    check(speedup > 5.0, f"expected a large speedup, got {speedup:.1f}x")
    print(f"    L={L} T={T} band w=100 (span {band.max_width}): "
          f"numpy {t1-t0:.3f}s, numba {t2-t1:.4f}s -> {speedup:.1f}x")
    per_locus = (t2 - t1) * 500 * 2 / 60
    print(f"    extrapolated: 500 reads x 2 assemblies = {per_locus:.2f} min/locus")


def main():
    tests = [
        test_numba_is_available,
        test_matches_the_oracle,
        test_terminal_row_matches,
        test_banded_matches,
        test_stable_on_a_long_read,
        test_speedup,
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
