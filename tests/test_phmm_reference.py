"""Tests for the reference forward and Viterbi (Stage B of docs/phmm.md).

The anchor test is `test_exact_path_enumeration`: for tiny inputs it walks the
transition graph of model.tex §5.3 exhaustively and sums whole paths, then
demands the DP match to 1e-12. Everything else in the pair HMM is that
recursion applied at scale, so ground truth here is worth more than any number
of plausibility checks.

Run: venv/bin/python tests/test_phmm_reference.py
"""

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qamodel.phmm import PhmmParams, quality_table
from qamodel.phmm.reference import forward, forward_window, viterbi

NEG_INF = -math.inf
QT = quality_table()

_failures = []


def check(cond, msg):
    if not cond:
        _failures.append(msg)
        print(f"  FAIL: {msg}")
    return cond


# --------------------------------------------------------------------------
# Ground truth: exhaustive walk of the transition graph.
# --------------------------------------------------------------------------

def enumerate_paths(read, qual, template, params, qt, origin=None,
                    full_template=False):
    """Log-sum over every legal state path, by brute force.

    Independent of reference.py by construction: it walks the §5.3 edge list
    directly rather than reusing any recursion.
    """
    tr = params.log_trans()
    L, T = len(read), len(template)
    read_u, tmpl_u = read.upper(), template.upper()
    out = []

    def emit_m(u, v):
        q = int(qual[u])
        return qt.log_match[q] if read_u[u] == tmpl_u[v] else qt.log_mismatch[q]

    def walk(u, v, state, logp):
        """u read bases emitted, v template consumed, currently in `state`."""
        if logp == NEG_INF:
            return
        if u == L:
            # D -> E is forbidden, so only M and I may finish a path.
            if state in ("M", "I") and (not full_template or v == T):
                out.append(logp)
            return
        if state == "M":
            opts = (("M", tr.m_m), ("I", tr.m_i), ("D", tr.m_d))
        elif state == "I":
            opts = (("I", tr.i_i), ("M", tr.i_m))
        else:
            opts = (("D", tr.d_d), ("M", tr.d_m))
        for nxt, lt in opts:
            if lt == NEG_INF:
                continue
            if nxt == "M":
                if v < T:
                    walk(u + 1, v + 1, "M", logp + lt + emit_m(u, v))
            elif nxt == "I":
                walk(u + 1, v, "I", logp + lt + qt.log_insert)
            else:
                if v < T:
                    walk(u, v + 1, "D", logp + lt)

    for o in (range(T) if origin is None else [origin]):
        # B -> M consumes t[o]; B -> I consumes nothing, so v stays at o.
        walk(1, o + 1, "M", tr.b_m + emit_m(0, o))
        if tr.b_i != NEG_INF:
            walk(1, o, "I", tr.b_i + qt.log_insert)

    if not out:
        return NEG_INF
    m = max(out)
    return m + math.log(sum(math.exp(x - m) for x in out))


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


# --------------------------------------------------------------------------

def test_exact_path_enumeration():
    """The anchor: DP must equal an exhaustive sum over paths."""
    rng = random.Random(20260817)
    params = PhmmParams(delta_i=0.02, delta_d=0.03, gamma_i=0.25, gamma_d=0.35)
    cases = 0
    worst = 0.0
    for _ in range(12):
        L = rng.randint(1, 4)
        T = rng.randint(1, 6)
        read, tmpl = _rand_seq(L, rng), _rand_seq(T, rng)
        qual = [rng.choice((20, 30, 40)) for _ in range(L)]
        for origin in (None, 0, min(1, T - 1)):
            for full in (False, True):
                want = enumerate_paths(read, qual, tmpl, params, QT, origin, full)
                got = forward(read, qual, tmpl, params, QT,
                              origin=origin, full_template=full).logp
                cases += 1
                if want == NEG_INF or got == NEG_INF:
                    check(want == got,
                          f"impossible case disagrees: {want} vs {got} "
                          f"({read}/{tmpl} origin={origin} full={full})")
                    continue
                worst = max(worst, abs(want - got))
                check(abs(want - got) < 1e-12,
                      f"forward {got!r} != enumerated {want!r} "
                      f"(read={read} tmpl={tmpl} origin={origin} full={full})")
    print(f"    {cases} configurations, worst |delta| = {worst:.2e}")


def test_window_additivity():
    """Free start must equal the log-sum of the per-origin passes.

    This is the definition of correctness for the windowing, and the property
    that makes an origin-anchored kernel a legitimate substitute at stage F.
    """
    rng = random.Random(11)
    params = PhmmParams(delta_i=0.01, delta_d=0.01)
    read = _rand_seq(12, rng)
    tmpl = _rand_seq(30, rng)
    qual = [35] * len(read)

    per_origin = [forward(read, qual, tmpl, params, QT, origin=o).logp
                  for o in range(len(tmpl))]
    finite = [x for x in per_origin if math.isfinite(x)]
    m = max(finite)
    want = m + math.log(sum(math.exp(x - m) for x in finite))
    got = forward_window(read, qual, tmpl, params, QT).logp
    check(abs(want - got) < 1e-9,
          f"free start {got!r} != logsumexp over origins {want!r}")


def test_homopolymer_degeneracy():
    """model.tex Remark 6, corrected: the ratio is 5, not 7.

    Template AAAAAAA against read AAAAAA under global alignment. One base is
    deleted and seven template positions could host that deletion, but position
    1 needs B -> D and position 7 needs D -> E, both forbidden. The five
    survivors each carry four M->M transitions regardless of where the deletion
    sits, so they are equiprobable and the forward/Viterbi ratio is exact.

    delta_i = 0 removes insertion paths, which would otherwise contribute
    O(delta_i * delta_d) and make the ratio 5 only to a few decimals.
    """
    params = PhmmParams(delta_i=0.0, delta_d=1e-3, gamma_d=0.3)
    kw = dict(params=params, qt=QT, origin=0, full_template=True)
    f = forward("AAAAAA", 40, "AAAAAAA", **kw).logp
    v = viterbi("AAAAAA", 40, "AAAAAAA", **kw).logp
    ratio = math.exp(f - v)
    check(abs(ratio - 5.0) < 1e-9,
          f"forward/Viterbi must be exactly 5, got {ratio!r}")

    vr = viterbi("AAAAAA", 40, "AAAAAAA", **kw)
    check(vr.path.count("D") == 1 and vr.path.count("M") == 6,
          f"best path must be 6 M and 1 D, got {vr.path}")
    check(not vr.path.startswith("D") and not vr.path.endswith("D"),
          f"forbidden transitions used: {vr.path}")


def test_mismatch_costs_s():
    """One substitution must cost exactly s = log(3(1-eps)/eps)."""
    params = PhmmParams(delta_i=0.0, delta_d=0.0)
    tmpl = "ACGTACGTAC"
    kw = dict(params=params, qt=QT, origin=0, full_template=True)
    perfect = forward(tmpl, 30, tmpl, **kw).logp
    mutated = "ACGTTCGTAC"          # position 4 A -> T
    one_off = forward(mutated, 30, tmpl, **kw).logp
    s = float(QT.cost()[30])
    check(abs((perfect - one_off) - s) < 1e-9,
          f"one mismatch cost {perfect - one_off:.6f}, expected s = {s:.6f}")


def test_monotone_in_edit_distance():
    params = PhmmParams(delta_i=0.0, delta_d=0.0)
    tmpl = "ACGTACGTACGTACGT"
    kw = dict(params=params, qt=QT, origin=0, full_template=True)
    scores = []
    read = list(tmpl)
    for i in (0, 3, 6, 9):
        read[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[read[i]]
        scores.append(forward("".join(read), 30, tmpl, **kw).logp)
    check(all(a > b for a, b in zip(scores, scores[1:])),
          f"log-likelihood must fall with each mismatch: {scores}")


def test_copy_number_signal():
    """Three exact tandem copies must score log 3 above one.

    Proposition 4's mechanism in miniature. This is a standing regression test:
    if it breaks, the method is blind to collapsed arrays.
    """
    rng = random.Random(7)
    unit = _rand_seq(40, rng)
    # Gaps off, so the only paths are the three exact placements and the score
    # is exactly log 3 rather than log 3 + O(delta).
    params = PhmmParams(delta_i=0.0, delta_d=0.0)
    one = forward_window(unit, 40, unit, params, QT).logp
    three = forward_window(unit, 40, unit * 3, params, QT).logp
    check(abs((three - one) - math.log(3.0)) < 1e-12,
          f"3 copies scored {three - one:.12f} above 1, expected log 3 = "
          f"{math.log(3.0):.12f}")


def test_padding_is_exact_without_gaps():
    """With gaps disabled, flanking template is unreachable and the score is
    invariant to machine precision.

    This is the clean statement of what the two exclusions buy. Each origin is
    counted exactly once, so padding neither duplicates a placement nor invents
    one; the added origins can only explain a random read through ~L mismatches,
    which is exp(-200) and contributes nothing measurable.
    """
    rng = random.Random(3)
    read = _rand_seq(25, rng)
    params = PhmmParams(delta_i=0.0, delta_d=0.0)
    base = forward_window(read, 40, read, params, QT).logp
    for pad in (1, 5, 20, 60):
        tmpl = _rand_seq(pad, rng) + read + _rand_seq(pad, rng)
        got = forward_window(read, 40, tmpl, params, QT).logp
        check(abs(got - base) < 1e-12,
              f"pad={pad} changed the gapless score by {got - base:.3e}")


def test_padding_adds_only_bounded_gap_mass():
    """With gaps on, padding does move the score, and it should.

    A pad that extends a homopolymer across the boundary lets the read's last
    base match one, two or three positions further right after deletions --
    Remark 6's degeneracy in a second guise, worth (1 - gamma_D^k) * delta_D.
    The increment is therefore data-dependent, but it is always non-negative
    (padding only ever adds paths) and bounded by the total gap-open mass.
    """
    params = PhmmParams()
    bound = (params.delta_d / (1.0 - params.gamma_d)
             + params.delta_i / (1.0 - params.gamma_i))
    rng = random.Random(3)
    read = _rand_seq(25, rng)
    base = forward_window(read, 40, read, params, QT).logp

    cases = [("homopolymer", read + "T" * 8),
             ("random right", read + _rand_seq(40, rng)),
             ("random left", _rand_seq(40, rng) + read)]
    for label, tmpl in cases:
        delta = forward_window(read, 40, tmpl, params, QT).logp - base
        check(delta >= -1e-12, f"{label}: padding removed mass ({delta:.3e})")
        check(delta <= bound,
              f"{label}: padding added {delta:.3e}, above the gap-open "
              f"bound {bound:.3e}")


def test_termination_is_a_uniform_prefactor():
    """tau > 0 must shift the score by exactly log tau, independent of length.

    Verifies eq. (9) carries termination out front, and that no path terminates
    twice. The residual tau inside M->M and I->M is removed here by holding the
    gap parameters at zero, isolating the prefactor.
    """
    rng = random.Random(5)
    tmpl = _rand_seq(35, rng)
    for L in (8, 20, 35):
        read = tmpl[:L]
        a = forward(read, 40, tmpl, PhmmParams(delta_i=0.0, delta_d=0.0, tau=0.0),
                    QT, origin=0).logp
        b = forward(read, 40, tmpl, PhmmParams(delta_i=0.0, delta_d=0.0, tau=1e-3),
                    QT, origin=0).logp
        # With no gaps the only tau dependence left is M->M = 1 - tau, applied
        # L-1 times, plus the prefactor the conditioned form omits.
        expected = (L - 1) * math.log1p(-1e-3)
        check(abs((b - a) - expected) < 1e-9,
              f"L={L}: tau shifted the score by {b - a:.9f}, expected {expected:.9f}")


def test_viterbi_bounded_by_forward():
    rng = random.Random(99)
    params = PhmmParams(delta_i=0.01, delta_d=0.01)
    for _ in range(6):
        read, tmpl = _rand_seq(rng.randint(4, 12), rng), _rand_seq(rng.randint(6, 20), rng)
        f = forward_window(read, 33, tmpl, params, QT).logp
        v = viterbi(read, 33, tmpl, params, QT).logp
        check(v <= f + 1e-12, f"Viterbi {v} exceeded forward {f}")


def test_viterbi_locates_the_read():
    """Origin and end must point at where the read actually sits."""
    rng = random.Random(42)
    read = _rand_seq(30, rng)
    left = _rand_seq(17, rng)
    tmpl = left + read + _rand_seq(23, rng)
    vr = viterbi(read, 40, tmpl, PhmmParams(), QT)
    check(vr.origin == 17, f"origin should be 17, got {vr.origin}")
    check(vr.end == 17 + 30 - 1, f"end should be 46, got {vr.end}")
    check(vr.path == "M" * 30, f"expected 30 matches, got {vr.path}")


def test_edge_mass_flags_a_truncated_window():
    """A read running off the right edge must be reported, not scored silently."""
    rng = random.Random(64)
    read = _rand_seq(40, rng)
    inside = forward_window(read, 40, _rand_seq(15, rng) + read + _rand_seq(15, rng),
                            PhmmParams(), QT)
    truncated = forward_window(read, 40, _rand_seq(15, rng) + read[:25],
                               PhmmParams(), QT)
    check(inside.edge_mass(margin=2) < 1e-6,
          f"a contained read must show no edge mass, got {inside.edge_mass(margin=2):.3e}")
    check(truncated.edge_mass(margin=2) > 0.5,
          f"a read hanging off the edge must show edge mass, got "
          f"{truncated.edge_mass(margin=2):.3e}")


def main():
    tests = [
        test_exact_path_enumeration,
        test_window_additivity,
        test_homopolymer_degeneracy,
        test_mismatch_costs_s,
        test_monotone_in_edit_distance,
        test_copy_number_signal,
        test_padding_is_exact_without_gaps,
        test_padding_adds_only_bounded_gap_mass,
        test_termination_is_a_uniform_prefactor,
        test_viterbi_bounded_by_forward,
        test_viterbi_locates_the_read,
        test_edge_mass_flags_a_truncated_window,
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
