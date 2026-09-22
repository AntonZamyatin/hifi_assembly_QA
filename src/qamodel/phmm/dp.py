"""The pair HMM dynamic program (model.tex §5.5).

Two engines compute the same recursion and must return the same number:

  _dp_full_table   every cell of (L+1) x (T+1), both semirings, back-pointers.
                   Written to be obviously the recursion rather than to be fast,
                   and it keeps the scalar deletion loop for that reason. It is
                   the oracle, and it is what `forward` and `viterbi` use.

  _dp_two_rows     two rows instead of the whole table, a vectorised deletion
                   row, and an optional band. Forward only -- traceback needs
                   every row's pointers. It is what `forward_window` uses, and
                   it is the only one that can run on a real read.

Public surface, three functions:

  forward(...)         oracle; carries `full_template` for the tests
  viterbi(...)         best single path, with traceback; never used to score
  forward_window(...)  production scoring, row-recycled and bandable

Coordinates.  The read is `r[0..L-1]` and the template `t[0..T-1]`.  The DP runs
over `(u, v)` with `u` read bases emitted and `v` the *absolute* template index
reached, so `v` is the number of template bases consumed counting from position
zero of the window, not from the read's origin.  Using absolute indices is what
lets a single matrix hold every origin at once.

  f_M[u, v]  last state M, aligning r[u-1] to t[v-1]
  f_I[u, v]  last state I, having emitted r[u-1] and consumed t[0..v-1]
  f_D[u, v]  last state D, consuming t[v-1]

Two extents, not one.  A window has a *template* extent, which must be long
enough to hold the read plus indel slack, and an *origin* extent, the set of
start positions the window is responsible for.  They are not the same interval
and conflating them is a correctness bug, so `origin` controls the second
independently:

  origin=None      every position in the template  (whole-unit free-start pass)
  origin=(lo, hi)  half-open range of origins, template may extend past `hi`
  origin=o         exactly one origin

The model's sum runs over every origin in the recruited *unit*, and a read in a
repeat has several placements, so a single call is in general one term of that
sum rather than the whole of it.  Totals are combined by log-sum-exp over
windows, which is valid only when their origin ranges are disjoint: template
extents overlap by design, origin extents must not, or shared origins are
counted twice.  `origin=None` over the whole unit is the exact case where one
call is the entire sum, and it is what stage F is checked against.

Ends.  `full_template=False` is free-end: the score sums terminal mass over
every template consumption, which absorbs trailing padding.  `f_D` is excluded
from that sum, and that exclusion *is* `D -> E` being forbidden.  There is no
analogous "end extent" to bound, because a path that runs off the right of the
window simply carries no mass there; `WindowScore.edge_mass` reports when that
has happened.

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


def _qual_array(qual, length: int) -> np.ndarray:
    """Accept a per-base quality array or one value broadcast over the read."""
    arr = np.full(length, int(qual), dtype=np.int64) if np.isscalar(qual) \
        else np.asarray(qual, dtype=np.int64)
    if arr.shape != (length,):
        raise ValueError(f"qual has length {arr.shape}, read has length {length}")
    return arr


def _origin_range(origin, T: int) -> tuple[int, int]:
    """Normalise `origin` to a half-open range of legal start positions.

    Origin `T` is never legal: it would put the read's start past the last
    template base.
    """
    if origin is None:
        return 0, T
    lo, hi = (origin, origin + 1) if isinstance(origin, (int, np.integer)) \
        else (int(origin[0]), int(origin[1]))
    if not 0 <= lo < hi <= T:
        raise ValueError(
            f"origin range [{lo}, {hi}) outside template of length {T}")
    return lo, hi


def _log_add(a: float, b: float) -> float:
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
    band_edge: float = 0.0

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


def _dp_full_table(read: np.ndarray, qual: np.ndarray, tmpl: np.ndarray,
        params: PhmmParams, qt: QualityTable,
        origin, use_max: bool, want_ptr: bool):
    """Shared recursion. `use_max` switches the semiring from sum to max.

    Forward and Viterbi share one body deliberately: the homopolymer test
    compares them, so a divergence between two separate implementations would
    show up as a modelling result rather than as a bug.
    """
    tr = params.log_trans()
    L, T = read.shape[0], tmpl.shape[0]
    lo, hi = _origin_range(origin, T)

    combine = max if use_max else _log_add

    def combine_vec(*arrays):
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

    # Begin state, seeded across the origin range only. B -> M at origin o aligns
    # r[0] to t[o] and so lands at v = o+1; B -> I emits r[0] without consuming
    # template and so stays at v = o. That one-position offset is why the two
    # slices differ.
    start_m = np.full(T + 1, NEG_INF)
    start_i = np.full(T + 1, NEG_INF)
    start_m[lo + 1:hi + 1] = tr.b_m
    start_i[lo:hi] = tr.b_i

    for u in range(1, L + 1):
        q = int(qual[u - 1])
        base = read[u - 1]
        hit = (tmpl == base) & (base != _N_CODE) & (tmpl != _N_CODE)
        emit_m = np.where(hit, qt.log_match[q], qt.log_mismatch[q])

        # M and I depend only on row u-1, so both vectorise over v.
        cand_m = (fM[u - 1, :T] + tr.m_m, fI[u - 1, :T] + tr.i_m, fD[u - 1, :T] + tr.d_m)
        inc_m = combine_vec(*cand_m)
        if u == 1:
            inc_m = combine_vec(inc_m, start_m[1:])
        fM[u, 1:] = emit_m + inc_m

        cand_i = (fM[u - 1] + tr.m_i, fI[u - 1] + tr.i_i)
        inc_i = combine_vec(*cand_i)
        if u == 1:
            inc_i = combine_vec(inc_i, start_i)
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


def _deletion_row(cur_m: np.ndarray, m_d: float, d_d: float) -> np.ndarray:
    """One D row, without the left-to-right scalar scan.

    The recurrence is `D[v] = LSE( M[v-1] + m_d, D[v-1] + d_d )`. Unrolled, it
    says a deletion reaching `v` opened at some `k < v` and was extended `v-1-k`
    times:

        D[v] = LSE over k < v of ( M[k] + m_d + (v-1-k)*d_d )

    The `(v-1)*d_d` part depends only on `v`, so pull it out; what is left is a
    running total over `k`, which is one `logaddexp.accumulate` instead of `T`
    scalar steps. The intermediate terms reach `T*|d_d|` in magnitude -- about
    18000 at `T = 15000`, `gamma_D = 0.3` -- so the reconstruction loses roughly
    `T*|d_d|*2^-52` of absolute accuracy, around 4e-12 at that size. Fine for the
    1e-9 tolerances used here, but it is why `_dp_full_table` keeps the scalar loop
    and why this is tested against it rather than trusted.
    """
    n = cur_m.shape[0]
    out = np.full(n, NEG_INF)
    if m_d == NEG_INF or n < 2:
        return out
    if d_d == NEG_INF:
        out[1:] = cur_m[:-1] + m_d          # gaps cannot extend: length 1 only
        return out
    k = np.arange(n, dtype=np.float64)
    c = np.logaddexp.accumulate(cur_m - k * d_d)
    out[1:] = m_d + (k[1:] - 1.0) * d_d + c[:-1]
    return out


def _logsumexp(*arrays) -> float:
    """log-sum-exp over several small arrays, ignoring -inf."""
    vals = np.concatenate([a[np.isfinite(a)] for a in arrays]) if arrays else None
    if vals is None or vals.size == 0:
        return NEG_INF
    m = vals.max()
    return float(m + np.log(np.exp(vals - m).sum()))


def _gather(src, src_lo, src_hi, want_lo, n, out):
    """`out[j] = src[want_lo + j]`, or -inf where that is outside the source band.

    Rows are held in absolute template coordinates, so cells outside a row's own
    band hold stale values from two rows ago. Clamping reads to the source band
    is what keeps that staleness unreachable, and it is cheaper than clearing.
    """
    out[:] = NEG_INF
    a = max(want_lo, src_lo)
    b = min(want_lo + n, src_hi)
    if b > a:
        out[a - want_lo:b - want_lo] = src[a:b]
    return out


def _dp_two_rows(read: np.ndarray, qual: np.ndarray, tmpl: np.ndarray,
                  params: PhmmParams, qt: QualityTable,
                  origin, band=None):
    """Forward pass keeping two rows, optionally restricted to a band.

    Row `u` reads only row `u-1`, so once row `u` exists row `u-2` is dead. Two
    buffers therefore suffice, swapped each step: memory drops from `O(L*T)` to
    `O(T)`, which is 9 KB rather than 5.1 GB for a 15 kb read against a 15.1 kb
    template. Scoring needs only the final row, so nothing is lost.

    Viterbi cannot do this: traceback walks backwards and needs every row's
    pointers, which is why `_dp_full_table` still allocates the whole thing.

    With a band, only `lo[u] <= v < hi[u]` is evaluated. Paths leaving the band
    are lost rather than approximated, which is why the second return value
    reports the largest share of any row's mass found sitting on a band edge.
    Edges that coincide with the template's own ends do not count: mass there is
    the read genuinely reaching the boundary, not the band cutting it off.

    Returns `(terminal_row, band_edge_fraction)`.
    """
    tr = params.log_trans()
    L, T = read.shape[0], tmpl.shape[0]
    lo_o, hi_o = _origin_range(origin, T)

    if band is None:
        band_lo = np.zeros(L + 1, dtype=np.int64)
        band_hi = np.full(L + 1, T + 1, dtype=np.int64)
    else:
        if band.rows != L + 1:
            raise ValueError(
                f"band has {band.rows} rows, read of length {L} needs {L + 1}")
        band_lo, band_hi = band.lo, band.hi

    start_m = np.full(T + 1, NEG_INF)
    start_i = np.full(T + 1, NEG_INF)
    start_m[lo_o + 1:hi_o + 1] = tr.b_m
    start_i[lo_o:hi_o] = tr.b_i

    prev_m = np.full(T + 1, NEG_INF)
    prev_i = np.full(T + 1, NEG_INF)
    prev_d = np.full(T + 1, NEG_INF)
    cur_m = np.full(T + 1, NEG_INF)
    cur_i = np.full(T + 1, NEG_INF)
    cur_d = np.full(T + 1, NEG_INF)

    w = int((band_hi - band_lo).max())
    buf = [np.empty(w) for _ in range(5)]

    p_lo, p_hi = int(band_lo[0]), int(band_hi[0])
    band_edge = 0.0

    for u in range(1, L + 1):
        a, b = int(band_lo[u]), int(band_hi[u])
        n = b - a
        am = max(a, 1)                      # M has no cell at v = 0
        nm = b - am

        q = int(qual[u - 1])
        base = read[u - 1]

        gm = _gather(prev_m, p_lo, p_hi, am - 1, nm, buf[0][:nm])
        gi = _gather(prev_i, p_lo, p_hi, am - 1, nm, buf[1][:nm])
        gd = _gather(prev_d, p_lo, p_hi, am - 1, nm, buf[2][:nm])
        inc_m = np.logaddexp(np.logaddexp(gm + tr.m_m, gi + tr.i_m), gd + tr.d_m)

        hm = _gather(prev_m, p_lo, p_hi, a, n, buf[3][:n])
        hi_ = _gather(prev_i, p_lo, p_hi, a, n, buf[4][:n])
        inc_i = np.logaddexp(hm + tr.m_i, hi_ + tr.i_i)

        if u == 1:
            inc_m = np.logaddexp(inc_m, start_m[am:b])
            inc_i = np.logaddexp(inc_i, start_i[a:b])

        tsl = tmpl[am - 1:b - 1]
        hit = (tsl == base) & (base != _N_CODE) & (tsl != _N_CODE)
        cur_m[a:am] = NEG_INF
        cur_m[am:b] = np.where(hit, qt.log_match[q], qt.log_mismatch[q]) + inc_m
        cur_i[a:b] = qt.log_insert + inc_i
        cur_d[a:b] = _deletion_row(cur_m[a:b], tr.m_d, tr.d_d)

        if band is not None and n > 2:
            total = _logsumexp(cur_m[a:b], cur_i[a:b], cur_d[a:b])
            if math.isfinite(total):
                sides = []
                if a > 0:
                    sides += [cur_m[a:a + 1], cur_i[a:a + 1], cur_d[a:a + 1]]
                if b < T + 1:
                    sides += [cur_m[b - 1:b], cur_i[b - 1:b], cur_d[b - 1:b]]
                if sides:
                    band_edge = max(band_edge,
                                    float(np.exp(_logsumexp(*sides) - total)))

        prev_m, cur_m = cur_m, prev_m
        prev_i, cur_i = cur_i, prev_i
        prev_d, cur_d = cur_d, prev_d
        p_lo, p_hi = a, b

    terminal = np.full(T + 1, NEG_INF)
    terminal[p_lo:p_hi] = np.logaddexp(prev_m[p_lo:p_hi], prev_i[p_lo:p_hi])
    return terminal, band_edge


def _terminal_row(fM, fI, L: int, use_max: bool) -> np.ndarray:
    """Per-end-position mass, excluding f_D: that exclusion is `D -> E` forbidden."""
    return np.maximum(fM[L], fI[L]) if use_max else np.logaddexp(fM[L], fI[L])


def _restrict_to_end(terminal: np.ndarray, T: int) -> np.ndarray:
    """Keep only paths that finished at the template's right edge.

    This is `full_template`: classical global alignment when combined with
    `origin=0`. It belongs to the caller, not to the recursion -- the engines
    never needed to know about it.
    """
    out = np.full_like(terminal, NEG_INF)
    out[T] = terminal[T]
    return out


def _score_from_terminal(terminal: np.ndarray, use_max: bool) -> float:
    """Collapse the per-end-position row to one number."""
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
    return r, _qual_array(qual, r.shape[0]), t, params or PhmmParams(), qt or quality_table()


def forward(read, qual, template, params: PhmmParams | None = None,
            qt: QualityTable | None = None, *,
            origin=None, full_template: bool = False) -> WindowScore:
    """Sum over every alignment path, and over the origins in `origin`.

    `origin=None` covers the whole template, `(lo, hi)` a half-open range, an
    integer a single position.  Use the range form whenever the template is
    padded past the origins this call is responsible for -- summing several such
    calls is only valid if their ranges are disjoint.

    `full_template=True` additionally forces the path to finish at the right
    edge, giving classical global alignment when combined with `origin=0`.  That
    combination is for the oracle tests; scoring never uses it.
    """
    r, q, t, p, table = _prepare(read, qual, template, params, qt)
    fM, fI, _, _ = _dp_full_table(r, q, t, p, table, origin, False, False)
    terminal = _terminal_row(fM, fI, r.shape[0], False)
    if full_template:
        terminal = _restrict_to_end(terminal, t.shape[0])
    return WindowScore(logp=_score_from_terminal(terminal, False),
                       terminal=terminal)


def viterbi(read, qual, template, params: PhmmParams | None = None,
            qt: QualityTable | None = None, *,
            origin=None, full_template: bool = False) -> ViterbiResult:
    """Best single path, with traceback. Never used for scoring (Remark 6)."""
    r, q, t, p, table = _prepare(read, qual, template, params, qt)
    L, T = r.shape[0], t.shape[0]
    fM, fI, fD, (pM, pI, pD) = _dp_full_table(r, q, t, p, table, origin, True, True)

    terminal = _terminal_row(fM, fI, L, True)
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
                   qt: QualityTable | None = None, *, origin=None,
                   band=None) -> WindowScore:
    """Score a read against one window. This is the entry point for real work.

    Free start over `origin`, free end, no termination factor -- so with
    `origin=None` over a whole recruited unit it is the model's full origin sum,
    and over a padded window around one alignment it is a single term of that
    sum, to be combined with others by log-sum-exp over *disjoint* origin ranges.

    Same answer as `forward`, but in `O(T)` memory rather than `O(L*T)`, and
    optionally banded. `forward` allocates 5.1 GB for a 15 kb read against a
    15.1 kb template and simply cannot run; this needs about 9 KB.

    There is deliberately no `full_template` here. It forces paths to end at the
    template's right edge, which is right for the oracle tests and silently wrong
    for scoring -- it would discard every path that finishes before the template
    does, which is most of them. `forward` keeps it; real work does not want it.
    """
    r, q, t, p, table = _prepare(read, qual, template, params, qt)
    terminal, band_edge = _dp_two_rows(r, q, t, p, table, origin, band)
    return WindowScore(logp=_score_from_terminal(terminal, False),
                       terminal=terminal, band_edge=band_edge)
