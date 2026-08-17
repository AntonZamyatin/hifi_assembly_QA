"""Reference forward and Viterbi for the pair HMM (model.tex §5.5).

Stage B of docs/phmm.md.  Pure numpy, log space, unbanded, no optimisation.
This is the permanent oracle every later stage is checked against, so it is
written to be obviously the recursion of §5.5 rather than to be fast.

Coordinates.  The read is `r[0..L-1]` and the template `t[0..T-1]`.  The DP runs
over `(u, v)` with `u` read bases emitted and `v` the *absolute* template index
reached, so `v` is the number of template bases consumed counting from position
zero of the window, not from the read's origin.  Using absolute indices is what
lets a single matrix hold every origin at once.

  f_M[u, v]  last state M, aligning r[u-1] to t[v-1]
  f_I[u, v]  last state I, having emitted r[u-1] and consumed t[0..v-1]
  f_D[u, v]  last state D, consuming t[v-1]

Origins and ends.  With `origin=None` the pass is free-start: the begin state is
seeded at every template position, which is exactly the marginalisation over
origins of §7 and the reason `B -> D` must stay excluded (see params.py).  With
`full_template=False` it is free-end: the score sums terminal mass over every
template consumption, which absorbs trailing padding.  `f_D` is excluded from
that sum, and that exclusion *is* `D -> E` being forbidden.

Termination.  Under the default `tau = 0` no termination factor is applied: the
returned value is `log sum_m [f_M(L,m) + f_I(L,m)]`, the length-conditioned
quantity of §5.6.  A positive `tau` shifts every path by the same `log tau`, so
callers wanting the joint model over `(r, l)` add it themselves.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from qamodel.phmm.params import PhmmParams, QualityTable, quality_table

NEG_INF = -math.inf

_STATE_M, _STATE_I, _STATE_D, _STATE_START, _STATE_NONE = 0, 1, 2, 3, -1
_STATE_NAME = {_STATE_M: "M", _STATE_I: "I", _STATE_D: "D"}

_N_CODE = ord("N")


def encode(seq: str) -> np.ndarray:
    """Sequence to uint8 codes, upper-cased.

    Upper-casing is not cosmetic: reference FASTAs are routinely soft-masked, and
    comparing raw bytes would score every repeat-masked base as a mismatch.
    """
    return np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)


def _as_qual(qual, length: int) -> np.ndarray:
    """Accept a per-base quality array or one value broadcast over the read."""
    arr = np.full(length, int(qual), dtype=np.int64) if np.isscalar(qual) \
        else np.asarray(qual, dtype=np.int64)
    if arr.shape != (length,):
        raise ValueError(f"qual has length {arr.shape}, read has length {length}")
    return arr


def _la(a: float, b: float) -> float:
    """Scalar log-add-exp; `np.logaddexp` is correct here but ~5x slower."""
    if b > a:
        a, b = b, a
    if b == NEG_INF:
        return a
    return a + math.log1p(math.exp(b - a))


@dataclass(frozen=True)
class WindowScore:
    """Result of a forward pass.

    `terminal[v]` is the log mass of paths that emitted the whole read having
    consumed template through `t[v-1]`, so it is the distribution over *end*
    positions.  `logp` is its total (or its value at `v = T` under
    `full_template`).
    """

    logp: float
    terminal: np.ndarray

    def edge_mass(self, margin: int = 1) -> float:
        """Fraction of terminal mass ending within `margin` of the right edge.

        The trustworthiness check of docs/phmm.md: a non-negligible value means
        the read runs past the window and the origin sum was truncated.  This
        covers the right edge only; for the left, compare `viterbi().origin`
        against the same margin.  An exact origin distribution needs a backward
        pass, which stage F can add if the Viterbi proxy proves too coarse.
        """
        if not np.isfinite(self.logp):
            return 0.0
        lo = max(0, len(self.terminal) - margin)
        tail = self.terminal[lo:]
        finite = tail[np.isfinite(tail)]
        if finite.size == 0:
            return 0.0
        m = finite.max()
        return float(np.exp(m + np.log(np.exp(finite - m).sum()) - self.logp))


@dataclass(frozen=True)
class ViterbiResult:
    """Best single path: its score, where it started and ended, and its states.

    `origin` and `end` are 0-based inclusive template indices.  `end < origin`
    is possible only for a degenerate all-insertion path, which consumes no
    template at all.
    """

    logp: float
    origin: int
    end: int
    path: str


def _dp(read: np.ndarray, qual: np.ndarray, tmpl: np.ndarray,
        params: PhmmParams, qt: QualityTable,
        origin: int | None, full_template: bool,
        use_max: bool, want_ptr: bool):
    """Shared recursion. `use_max` switches the semiring from sum to max.

    Forward and Viterbi share one body deliberately: the homopolymer test
    compares them, so a divergence between two separate implementations would
    show up as a modelling result rather than as a bug.
    """
    tr = params.log_trans()
    L, T = read.shape[0], tmpl.shape[0]
    if origin is not None and not 0 <= origin < T:
        raise ValueError(f"origin {origin} outside template of length {T}")

    combine = max if use_max else _la

    def vcombine(*arrays):
        out = arrays[0]
        for a in arrays[1:]:
            out = np.maximum(out, a) if use_max else np.logaddexp(out, a)
        return out

    shape = (L + 1, T + 1)
    fM = np.full(shape, NEG_INF)
    fI = np.full(shape, NEG_INF)
    fD = np.full(shape, NEG_INF)

    pM = pI = pD = None
    if want_ptr:
        pM = np.full(shape, _STATE_NONE, dtype=np.int8)
        pI = np.full(shape, _STATE_NONE, dtype=np.int8)
        pD = np.full(shape, _STATE_NONE, dtype=np.int8)

    # Begin state. B -> M at origin o aligns r[0] to t[o], landing at v = o+1.
    # B -> I emits r[0] without consuming template, so it stays at v = o.
    start_m = np.full(T + 1, NEG_INF)
    start_i = np.full(T + 1, NEG_INF)
    if origin is None:
        start_m[1:] = tr.b_m
        start_i[:T] = tr.b_i
    else:
        start_m[origin + 1] = tr.b_m
        start_i[origin] = tr.b_i

    for u in range(1, L + 1):
        q = int(qual[u - 1])
        base = read[u - 1]
        hit = (tmpl == base) & (base != _N_CODE) & (tmpl != _N_CODE)
        emit_m = np.where(hit, qt.log_match[q], qt.log_mismatch[q])

        # M and I depend only on row u-1, so both vectorise over v.
        cand_m = (fM[u - 1, :T] + tr.m_m, fI[u - 1, :T] + tr.i_m, fD[u - 1, :T] + tr.d_m)
        inc_m = vcombine(*cand_m)
        if u == 1:
            inc_m = vcombine(inc_m, start_m[1:])
        fM[u, 1:] = emit_m + inc_m

        cand_i = (fM[u - 1] + tr.m_i, fI[u - 1] + tr.i_i)
        inc_i = vcombine(*cand_i)
        if u == 1:
            inc_i = vcombine(inc_i, start_i)
        fI[u] = qt.log_insert + inc_i

        if want_ptr:
            # Map argmax position -> state code explicitly; the two stacks have
            # different layouts and START sits at a different index in each.
            def _ptr(cands, extra, codes):
                stack = np.vstack(cands + extra)
                lookup = np.array(codes, dtype=np.int8)
                return np.where(np.isneginf(stack.max(axis=0)),
                                _STATE_NONE, lookup[stack.argmax(axis=0)])

            seeded = (u == 1)
            pM[u, 1:] = _ptr(cand_m, (start_m[1:],) if seeded else (),
                             [_STATE_M, _STATE_I, _STATE_D]
                             + ([_STATE_START] if seeded else []))
            pI[u] = _ptr(cand_i, (start_i,) if seeded else (),
                         [_STATE_M, _STATE_I]
                         + ([_STATE_START] if seeded else []))

        # D depends on the same row at v-1, so it is a left-to-right scan.
        row_m, row_d = fM[u], fD[u]
        prev = NEG_INF
        prev_ptr = _STATE_NONE
        for v in range(1, T + 1):
            from_m = row_m[v - 1] + tr.m_d
            from_d = prev + tr.d_d
            if want_ptr:
                prev_ptr = _STATE_M if from_m >= from_d else _STATE_D
                if from_m == NEG_INF and from_d == NEG_INF:
                    prev_ptr = _STATE_NONE
                pD[u, v] = prev_ptr
            prev = combine(from_m, from_d)
            row_d[v] = prev

    return fM, fI, fD, (pM, pI, pD)


def _terminal(fM, fI, L: int, use_max: bool) -> np.ndarray:
    """Per-end-position mass, excluding f_D: that exclusion is `D -> E` forbidden."""
    return np.maximum(fM[L], fI[L]) if use_max else np.logaddexp(fM[L], fI[L])


def _total(terminal: np.ndarray, full_template: bool, T: int, use_max: bool) -> float:
    if full_template:
        return float(terminal[T])
    finite = terminal[np.isfinite(terminal)]
    if finite.size == 0:
        return NEG_INF
    if use_max:
        return float(finite.max())
    m = finite.max()
    return float(m + np.log(np.exp(finite - m).sum()))


def _prepare(read, qual, template, params, qt):
    r = encode(read) if isinstance(read, str) else np.asarray(read, dtype=np.uint8)
    t = encode(template) if isinstance(template, str) else np.asarray(template, dtype=np.uint8)
    return r, _as_qual(qual, r.shape[0]), t, params or PhmmParams(), qt or quality_table()


def forward(read, qual, template, params: PhmmParams | None = None,
            qt: QualityTable | None = None, *,
            origin: int | None = None, full_template: bool = False) -> WindowScore:
    """Sum over every alignment path, and over origins unless one is given.

    `origin=None` marginalises over all template positions (the production
    case); an integer pins the start.  `full_template=True` additionally forces
    the path to finish at the right edge, giving classical global alignment when
    combined with `origin=0` -- used by the oracle tests, not in scoring.
    """
    r, q, t, p, table = _prepare(read, qual, template, params, qt)
    fM, fI, _, _ = _dp(r, q, t, p, table, origin, full_template, False, False)
    terminal = _terminal(fM, fI, r.shape[0], False)
    if full_template:
        masked = np.full_like(terminal, NEG_INF)
        masked[t.shape[0]] = terminal[t.shape[0]]
        terminal = masked
    return WindowScore(logp=_total(terminal, full_template, t.shape[0], False),
                       terminal=terminal)


def viterbi(read, qual, template, params: PhmmParams | None = None,
            qt: QualityTable | None = None, *,
            origin: int | None = None, full_template: bool = False) -> ViterbiResult:
    """Best single path, with traceback. Never used for scoring (Remark 6)."""
    r, q, t, p, table = _prepare(read, qual, template, params, qt)
    L, T = r.shape[0], t.shape[0]
    fM, fI, fD, (pM, pI, pD) = _dp(r, q, t, p, table, origin, full_template, True, True)

    terminal = _terminal(fM, fI, L, True)
    if full_template:
        best_v, best = T, float(terminal[T])
    else:
        finite = np.isfinite(terminal)
        if not finite.any():
            return ViterbiResult(NEG_INF, -1, -1, "")
        best_v = int(np.argmax(np.where(finite, terminal, NEG_INF)))
        best = float(terminal[best_v])
    if not math.isfinite(best):
        return ViterbiResult(NEG_INF, -1, -1, "")

    state = _STATE_M if fM[L, best_v] >= fI[L, best_v] else _STATE_I
    u, v, states = L, best_v, []
    end = v - 1
    while u > 0 or state == _STATE_D:
        states.append(_STATE_NAME[state])
        nxt = {_STATE_M: pM, _STATE_I: pI, _STATE_D: pD}[state][u, v]
        if state == _STATE_M:
            u, v = u - 1, v - 1
        elif state == _STATE_I:
            u = u - 1
        else:
            v = v - 1
        if nxt == _STATE_START or nxt == _STATE_NONE:
            break
        state = int(nxt)
    states.reverse()
    return ViterbiResult(logp=best, origin=v, end=end, path="".join(states))


def forward_window(read, qual, template, params: PhmmParams | None = None,
                   qt: QualityTable | None = None) -> WindowScore:
    """The production primitive: free start, free end, no termination factor."""
    return forward(read, qual, template, params, qt, origin=None, full_template=False)
