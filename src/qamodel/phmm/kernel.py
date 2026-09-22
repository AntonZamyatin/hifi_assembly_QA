"""Compiled dynamic program, third engine (stage F4 of docs/phmm.md).

Same recursion as `dp.py`, in linear space with per-row rescaling instead of log
space. That choice, not the compilation, is where most of the speed comes from.

Measured on the numpy engine at band width 200: `np.logaddexp` costs 3.66 us per
call against 0.41 us for `np.maximum`, and roughly 89% of that is the `exp` and
`log1p` pair. Across a full row those two operations are about 85% of the time.
In linear space they are gone entirely -- every cell is multiply-add:

    M = em * (mm*pM + im*pI + dm*pD)
    I = li * (mi*pM + ii*pI)
    D = md*curM + dd*prevD

about 13 flops per cell. The price is bookkeeping: each row is divided by its own
maximum and the log of that divisor is accumulated, so nothing overflows or
underflows however long the read. Cells more than ~700 nats below the row maximum
flush to zero, which is correct -- they contribute nothing.

This is why the earlier decision to stay in log space does not survive into a
compiled kernel. It was right for numpy, where the deletion scan had to be
vectorised through `logaddexp.accumulate` and a Python loop was unaffordable.
Inside numba a sequential loop is exactly what is wanted, and log space becomes
the wrong choice.

`dp.py` remains the oracle. This engine is checked against it, never trusted.
"""

from __future__ import annotations

import math

import numpy as np

from qamodel.phmm.dp import (NEG_INF, WindowScore, _origin_range, _prepare)
from qamodel.phmm.params import PhmmParams, QualityTable

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:                                     # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):                                  # type: ignore[misc]
        def wrap(fn):
            return fn
        return wrap if not a else a[0]


_N_CODE = 78            # ord('N'); numba wants a literal, not a module global lookup


@njit(cache=True)
def _run(read, qual, tmpl, band_lo, band_hi,
         lin_match, lin_mismatch, lin_insert,
         b_m, b_i, m_m, m_i, m_d, i_i, i_m, d_d, d_m,
         origin_lo, origin_hi):
    """Returns (row_m, row_i, log_scale, p_lo, p_hi, band_edge).

    `row_m + row_i` is the terminal row in the final row's own scale; the true
    log mass is `log(row_m[v] + row_i[v]) + log_scale`.
    """
    L = read.shape[0]
    T = tmpl.shape[0]
    n = T + 1

    pM = np.zeros(n); pI = np.zeros(n); pD = np.zeros(n)
    cM = np.zeros(n); cI = np.zeros(n); cD = np.zeros(n)

    log_scale = 0.0
    band_edge = 0.0
    p_lo = band_lo[0]
    p_hi = band_hi[0]

    for u in range(1, L + 1):
        a = band_lo[u]
        b = band_hi[u]
        q = qual[u - 1]
        base = read[u - 1]
        em_hit = lin_match[q]
        em_miss = lin_mismatch[q]
        seeded = (u == 1)

        # The buffer being written still holds row u-2 outside this row's band.
        for v in range(a, b):
            cM[v] = 0.0
            cI[v] = 0.0
            cD[v] = 0.0

        for v in range(a, b):
            vp = v - 1
            acc_m = 0.0
            if p_lo <= vp < p_hi:
                acc_m = m_m * pM[vp] + i_m * pI[vp] + d_m * pD[vp]
            if seeded and origin_lo <= vp < origin_hi:
                acc_m += b_m
            if v >= 1 and acc_m > 0.0:
                tb = tmpl[vp]
                if tb == base and base != _N_CODE and tb != _N_CODE:
                    cM[v] = em_hit * acc_m
                else:
                    cM[v] = em_miss * acc_m

            acc_i = 0.0
            if p_lo <= v < p_hi:
                acc_i = m_i * pM[v] + i_i * pI[v]
            if seeded and origin_lo <= v < origin_hi:
                acc_i += b_i
            cI[v] = lin_insert * acc_i

        # Deletions run left to right within the row: a gap is crossed without
        # advancing u, so this scan is inherently sequential. In linear space
        # each step is one multiply-add rather than an exp plus a log1p.
        run = 0.0
        for v in range(a, b):
            vp = v - 1
            mval = cM[vp] if vp >= a else 0.0
            run = m_d * mval + d_d * run
            cD[v] = run

        # Rescale by the row maximum so nothing drifts out of range.
        s = 0.0
        for v in range(a, b):
            if cM[v] > s:
                s = cM[v]
            if cI[v] > s:
                s = cI[v]
            if cD[v] > s:
                s = cD[v]
        if s <= 0.0:
            return cM, cI, -np.inf, a, b, band_edge
        inv = 1.0 / s
        total = 0.0
        for v in range(a, b):
            cM[v] *= inv
            cI[v] *= inv
            cD[v] *= inv
            total += cM[v] + cI[v] + cD[v]
        log_scale += math.log(s)

        if (b - a) > 2 and total > 0.0:
            edge = 0.0
            if a > 0:
                edge += cM[a] + cI[a] + cD[a]
            if b < n:
                edge += cM[b - 1] + cI[b - 1] + cD[b - 1]
            frac = edge / total
            if frac > band_edge:
                band_edge = frac

        tM = pM; pM = cM; cM = tM
        tI = pI; pI = cI; cI = tI
        tD = pD; pD = cD; cD = tD
        p_lo = a
        p_hi = b

    return pM, pI, log_scale, p_lo, p_hi, band_edge


def forward_numba(read, qual, template, params: PhmmParams | None = None,
                  qt: QualityTable | None = None, *,
                  origin=None, band=None) -> WindowScore:
    """`forward_window` on the compiled engine. Same arguments, same result."""
    if not HAVE_NUMBA:
        raise RuntimeError("numba is not installed; use forward_window instead")

    r, q, t, p, table = _prepare(read, qual, template, params, qt)
    L, T = r.shape[0], t.shape[0]
    lo_o, hi_o = _origin_range(origin, T)

    if band is None:
        band_lo = np.zeros(L + 1, dtype=np.int64)
        band_hi = np.full(L + 1, T + 1, dtype=np.int64)
    else:
        if band.rows != L + 1:
            raise ValueError(
                f"band has {band.rows} rows, read of length {L} needs {L + 1}")
        band_lo, band_hi = band.lo, band.hi

    tr = p.log_trans()
    row_m, row_i, log_scale, p_lo, p_hi, band_edge = _run(
        r, q, t, band_lo, band_hi,
        np.exp(table.log_match), np.exp(table.log_mismatch),
        math.exp(table.log_insert),
        math.exp(tr.b_m), math.exp(tr.b_i), math.exp(tr.m_m), math.exp(tr.m_i),
        math.exp(tr.m_d), math.exp(tr.i_i), math.exp(tr.i_m),
        math.exp(tr.d_d), math.exp(tr.d_m),
        lo_o, hi_o)

    terminal = np.full(T + 1, NEG_INF)
    if math.isfinite(log_scale):
        lin = row_m[p_lo:p_hi] + row_i[p_lo:p_hi]
        nz = lin > 0.0
        if nz.any():
            terminal[p_lo:p_hi][nz] = np.log(lin[nz]) + log_scale

    finite = terminal[np.isfinite(terminal)]
    if finite.size == 0:
        logp = NEG_INF
    else:
        mx = finite.max()
        logp = float(mx + np.log(np.exp(finite - mx).sum()))
    return WindowScore(logp=logp, terminal=terminal, band_edge=float(band_edge))
