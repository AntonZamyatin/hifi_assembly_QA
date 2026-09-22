"""Tests for placements -> windows (stage F3 of docs/phmm.md).

Everything here runs on synthetic sequence where the right answer is known, so
the strand and clip arithmetic is checked before any real read goes near it.
The yardstick throughout: a read that is an exact copy of some stretch of the
assembly must score, through its windows, the same as scoring it directly
against that stretch -- whichever strand it came from and however it is clipped.

Run: venv/bin/python tests/test_phmm_windows.py
"""

import math
import os
import random
import sys
from types import SimpleNamespace as NS

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qamodel.phmm import PhmmParams, quality_table
from qamodel.phmm.kernel import forward_numba
from qamodel.phmm.reads import revcomp
from qamodel.phmm.windows import Placement, build_windows, reverse_cigar

QT, P = quality_table(), PhmmParams()
_failures = []
rng = random.Random(20260918)


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


class FakeFasta:
    def __init__(self, contigs):
        self.c = contigs

    def fetch(self, contig, lo, hi):
        return self.c[contig][lo:hi]


def lse(xs):
    xs = [x for x in xs if math.isfinite(x)]
    if not xs:
        return -math.inf
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def via_windows(read, ws):
    return lse(forward_numba(read, 40, w.template, P, QT, origin=w.origin,
                             band=w.band).logp for w in ws.windows)


ASM = seq(6000)
FA = FakeFasta({"c1": ASM, "c2": seq(3000)})
SEG = [NS(contig="c1", start=500, end=5500)]
PERFECT = forward_numba(ASM[2000:3000], 40, ASM[2000:3000], P, QT, origin=0).logp


def test_reverse_cigar():
    check(reverse_cigar("5S10M2D3M1I") == "1I3M2D10M5S", reverse_cigar("5S10M2D3M1I"))


def test_forward_exact():
    read = ASM[2000:3000]
    ws = build_windows([Placement("c1", 2000, 3000, 1, "1000M", 0, 0)], SEG, 1000, FA)
    check(len(ws.windows) == 1 and not ws.straddles, f"{ws}")
    got = via_windows(read, ws)
    check(abs(got - PERFECT) < 1e-3, f"forward: {got:.5f} vs {PERFECT:.5f}")


def test_reverse_exact():
    """The read as sequenced is the reverse complement of the reference stretch."""
    read = revcomp(ASM[2000:3000])
    ws = build_windows([Placement("c1", 2000, 3000, 0, "1000M", 0, 0)], SEG, 1000, FA)
    got = via_windows(read, ws)
    check(abs(got - PERFECT) < 1e-3, f"reverse: {got:.5f} vs {PERFECT:.5f}")


def test_forward_soft_clip():
    """5 bp clipped at the read's start: the window must reach back for them."""
    read = ASM[1995:3000]
    ws = build_windows([Placement("c1", 2000, 3000, 1, "5S1000M", 5, 0)], SEG, 1005, FA)
    want = forward_numba(read, 40, read, P, QT, origin=0).logp
    got = via_windows(read, ws)
    check(abs(got - want) < 1e-3, f"forward clip: {got:.4f} vs {want:.4f}")


def test_reverse_soft_clip():
    """Clip recorded at the reference RIGHT end is the read's START when reversed.

    60 bp, deliberately larger than the k = 20 origin slack: a 5 bp clip put on
    the wrong side is silently absorbed by that slack, which is the design working
    but would make this test unable to tell right from wrong.
    """
    read = revcomp(ASM[2000:3060])
    ws = build_windows([Placement("c1", 2000, 3000, 0, "1000M60S", 0, 60)], SEG, 1060, FA)
    want = forward_numba(read, 40, read, P, QT, origin=0).logp
    got = via_windows(read, ws)
    check(abs(got - want) < 1e-3, f"reverse clip: {got:.4f} vs {want:.4f}")
    wrong = build_windows([Placement("c1", 2000, 3000, 0, "60S1000M", 60, 0)], SEG, 1060, FA)
    bad = via_windows(read, wrong)
    check(bad < want - 20.0,
          f"a 60 bp clip on the wrong side must cost, got {bad:.2f} vs {want:.2f}")


def test_straddle_and_outside():
    ws = build_windows([Placement("c1", 5400, 5600, 1, "200M", 0, 0)], SEG, 200, FA)
    check(ws.straddles and not ws.windows, "a placement crossing the segment end must straddle")
    ws = build_windows([Placement("c2", 100, 300, 1, "200M", 0, 0)], SEG, 200, FA)
    check(ws.outside == 1 and not ws.straddles and not ws.windows,
          "a placement on another contig is outside the unit")


def test_merge_keeps_origins_disjoint():
    near = [Placement("c1", 2000, 3000, 1, "1000M", 0, 0),
            Placement("c1", 2010, 3010, 1, "1000M", 0, 0)]
    ws = build_windows(near, SEG, 1000, FA, k=20)
    check(len(ws.windows) == 1 and ws.windows[0].members == 2,
          f"placements 10 bp apart must merge, got {len(ws.windows)} windows")
    far = [Placement("c1", 1000, 2000, 1, "1000M", 0, 0),
           Placement("c1", 3000, 4000, 1, "1000M", 0, 0)]
    ws = build_windows(far, SEG, 1000, FA, k=20)
    check(len(ws.windows) == 2, f"placements 2 kb apart stay separate, got {len(ws.windows)}")
    a, b = sorted(ws.windows, key=lambda w: w.ref_lo)
    ga = (a.ref_lo + a.origin[0], a.ref_lo + a.origin[1])
    gb = (b.ref_lo + b.origin[0], b.ref_lo + b.origin[1])
    check(ga[1] <= gb[0], f"owned starts overlap in reference coordinates: {ga} {gb}")


def test_tandem_copies_sum_to_log_3():
    """Three exact copies, three placements, three windows: log 3 above one copy."""
    unit = seq(800)
    asm = seq(500) + unit * 3 + seq(500)
    fa = FakeFasta({"t": asm})
    seg = [NS(contig="t", start=0, end=len(asm))]
    pl = [Placement("t", 500 + i * 800, 1300 + i * 800, 1, "800M", 0, 0) for i in range(3)]
    three = via_windows(unit, build_windows(pl, seg, 800, fa))
    one = via_windows(unit, build_windows(pl[:1], seg, 800, fa))
    check(abs((three - one) - math.log(3)) < 1e-6,
          f"three copies scored {three - one:.6f} above one, expected {math.log(3):.6f}")


def main():
    tests = [test_reverse_cigar, test_forward_exact, test_reverse_exact,
             test_forward_soft_clip, test_reverse_soft_clip,
             test_straddle_and_outside, test_merge_keeps_origins_disjoint,
             test_tandem_copies_sum_to_log_3]
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
