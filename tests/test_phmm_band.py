"""Tests for banding (F2 of docs/phmm.md).

Banding is a restriction, not an approximation with a knob: paths that leave the
band are *lost*, not down-weighted. So the tests come in pairs -- a generous band
must reproduce the unbanded answer exactly, and a tight one must both lose mass
and say so through `band_edge`. A band that silently loses mass is the failure
mode worth guarding, because it looks like a real likelihood difference.

Run: venv/bin/python tests/test_phmm_band.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from qamodel.phmm import PhmmParams, quality_table
from qamodel.phmm.band import Band, cigar_path
from qamodel.phmm.dp import forward_window

QT = quality_table()
_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_cigar_path_geometry():
    """v[u] is the template position reached after u read bases."""
    check(list(cigar_path(6, "50M", 10)) == [10, 11, 12, 13, 14, 15, 16],
          "plain matches must walk the diagonal from ref_start")
    check(list(cigar_path(6, "3M2D3M", 0)) == [0, 1, 2, 5, 6, 7, 8],
          "a deletion must advance v without advancing u")
    check(list(cigar_path(8, "3M2I3M", 0)) == [0, 1, 2, 3, 3, 3, 4, 5, 6],
          "an insertion must hold v still")
    # A leading clip backs the path off the reported start; otherwise the first
    # matched base lands `clip` positions too far right.
    check(list(cigar_path(10, "5S5M", 100))[5] == 100,
          f"5S5M at 100: read base 5 must align to template 100, "
          f"got {cigar_path(10, '5S5M', 100)[5]}")
    check(list(cigar_path(10, "5M5S", 100))[:6] == [100, 101, 102, 103, 104, 105],
          "a trailing clip must continue the diagonal")
    for bad in ("", "10Z", "M10", "10 M"):
        try:
            cigar_path(10, bad)
        except ValueError:
            continue
        _failures.append(f"CIGAR {bad!r} must be rejected")
        print(f"  FAIL: CIGAR {bad!r} must be rejected")


def test_band_rejects_empty_rows():
    try:
        Band(np.array([5, 5]), np.array([5, 6]))
    except ValueError:
        return
    _failures.append("an empty band row must raise")
    print("  FAIL: an empty band row must raise")


def test_generous_band_reproduces_unbanded():
    """A band wider than any drift must change nothing at all."""
    rng = random.Random(3)
    worst = 0.0
    for _ in range(8):
        T = rng.randint(40, 400)
        tmpl = _rand_seq(T, rng)
        start = rng.randint(0, T // 3)
        L = rng.randint(20, max(21, T // 2))
        read = tmpl[start:start + L]
        params = PhmmParams()
        ref = forward_window(read, 38, tmpl, params, QT,
                           origin=(start, start + 1)).logp
        band = Band.diagonal(len(read), T, start, start + 1, width=T)
        got = forward_window(read, 38, tmpl, params, QT,
                           origin=(start, start + 1), band=band)
        worst = max(worst, abs(ref - got.logp))
        check(abs(ref - got.logp) < 1e-9,
              f"generous band differs by {abs(ref-got.logp):.2e}")
        check(got.band_edge < 1e-6,
              f"generous band should report no edge mass, got {got.band_edge:.2e}")
    print(f"    worst |delta| = {worst:.2e}")


def test_tight_band_loses_mass_and_reports_it():
    """The guarded failure mode: a band too narrow must not fail silently."""
    rng = random.Random(21)
    tmpl = _rand_seq(600, rng)
    read = tmpl[100:400]
    read = read[:150] + _rand_seq(40, rng) + read[150:]   # 40 bp insertion
    params = PhmmParams(gamma_i=0.9)
    L = len(read)
    ref = forward_window(read, 38, tmpl, params, QT, origin=(100, 101)).logp

    wide = forward_window(read, 38, tmpl, params, QT, origin=(100, 101),
                        band=Band.diagonal(L, 600, 100, 101, width=100))
    tight = forward_window(read, 38, tmpl, params, QT, origin=(100, 101),
                         band=Band.diagonal(L, 600, 100, 101, width=5))
    check(abs(wide.logp - ref) < 1e-9,
          f"width 100 must cover a 40 bp indel; off by {abs(wide.logp-ref):.2e}")
    check(tight.logp < ref - 10.0,
          f"width 5 must lose real mass: {tight.logp:.2f} vs {ref:.2f}")
    check(tight.band_edge > 0.01,
          f"width 5 must report edge mass, got {tight.band_edge:.2e}")
    print(f"    unbanded {ref:.3f}  w=100 {wide.logp:.3f} (edge {wide.band_edge:.1e})"
          f"  w=5 {tight.logp:.3f} (edge {tight.band_edge:.2f})")


def test_cigar_band_follows_a_large_deletion():
    """Why the CIGAR option exists.

    A 300 bp deletion throws the path 300 columns off the straight diagonal, so
    a diagonal band would need width > 300 to contain it. A CIGAR-anchored band
    bends around the event and width 100 is plenty -- the width then only has to
    cover disagreement with the aligner's path, not the size of the event.
    """
    rng = random.Random(8)
    x, y, z = _rand_seq(500, rng), _rand_seq(300, rng), _rand_seq(500, rng)
    tmpl = x + y + z
    read = x + z
    params = PhmmParams(gamma_d=0.9)
    L, T = len(read), len(tmpl)
    ref = forward_window(read, 40, tmpl, params, QT, origin=(0, 1)).logp

    diag = forward_window(read, 40, tmpl, params, QT, origin=(0, 1),
                        band=Band.diagonal(L, T, 0, 1, width=100))
    cig = forward_window(read, 40, tmpl, params, QT, origin=(0, 1),
                       band=Band.from_cigar(L, T, "500M300D500M", 0, width=100))

    check(diag.logp < ref - 50.0,
          f"a width-100 diagonal must miss a 300 bp deletion; "
          f"got {diag.logp:.2f} vs {ref:.2f}")
    check(abs(cig.logp - ref) < 1e-9,
          f"the CIGAR band must recover it exactly; off by {abs(cig.logp-ref):.2e}")
    check(cig.band_edge < 1e-6,
          f"CIGAR band should sit clear of its edges, got {cig.band_edge:.2e}")
    print(f"    unbanded {ref:.3f}   diagonal w=100 {diag.logp:.3f}"
          f"   CIGAR w=100 {cig.logp:.3f}")


def test_banding_is_faster():
    rng = random.Random(11)
    T = 6000
    tmpl = _rand_seq(T, rng)
    read = tmpl[100:5900]
    L = len(read)
    band = Band.diagonal(L, T, 100, 101, width=100)

    t0 = time.perf_counter()
    slow = forward_window(read, 40, tmpl, PhmmParams(), QT, origin=(100, 101)).logp
    t1 = time.perf_counter()
    fast = forward_window(read, 40, tmpl, PhmmParams(), QT, origin=(100, 101),
                        band=band)
    t2 = time.perf_counter()
    check(abs(slow - fast.logp) < 1e-9,
          f"banded answer differs by {abs(slow-fast.logp):.2e}")
    speedup = (t1 - t0) / (t2 - t1)
    check(speedup > 3.0, f"expected a real speedup, got {speedup:.1f}x")
    print(f"    L={L} T={T}: unbanded {t1-t0:.2f}s, band w=100 {t2-t1:.2f}s "
          f"-> {speedup:.1f}x, width {band.max_width}")


def main():
    tests = [
        test_cigar_path_geometry,
        test_band_rejects_empty_rows,
        test_generous_band_reproduces_unbanded,
        test_tight_band_loses_mass_and_reports_it,
        test_cigar_band_follows_a_large_deletion,
        test_banding_is_faster,
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
