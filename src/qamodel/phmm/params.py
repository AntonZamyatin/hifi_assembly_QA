"""Pair HMM parameters and quality-derived emission constants.

Implements the transition table of model.tex §5.3 and the emissions of §5.4.
Stage A of docs/phmm.md.

Two structural facts about the model are carried by the types rather than by
comment.  `LogTrans` has no `d_e` field because termination from D is forbidden
(a trailing deletion consumes template past the last emitted base, which no
observation can distinguish from terminating earlier), and no `b_d` field
because B -> D duplicates the path from origin s+1 -- redundancy created by the
marginalisation over origins, not by the alignment model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields

import numpy as np

MAX_PHRED = 93          # SAM caps Phred at 93
EPS_MAX = 0.75          # at eps = 3/4 a mismatch is as likely as a match; beyond
                        # that the "error" model inverts and s goes negative
EPS_MIN = 1e-10         # guards log overflow at implausibly high reported quality


def _log(x: float) -> float:
    """log with -inf for zero, without numpy's divide-by-zero warning."""
    return math.log(x) if x > 0.0 else -math.inf


@dataclass(frozen=True)
class LogTrans:
    """Log transition probabilities, one field per edge of the §5.3 table."""

    b_m: float
    b_i: float
    m_m: float
    m_i: float
    m_d: float
    m_e: float
    i_i: float
    i_m: float
    i_e: float
    d_d: float
    d_m: float


@dataclass(frozen=True)
class PhmmParams:
    """Gap and termination probabilities, shared by both candidates.

    The defaults are placeholders of the right order for HiFi, to be replaced by
    the estimates of stage G.  They must never be fitted per candidate: doing so
    lets a poorer assembly compensate with looser error rates (model.tex §12.5).

    `tau` cancels from the comparison whenever the gap parameters are
    context-independent, because the normaliser Z is then origin- and
    candidate-independent (§5.6).  It matters only once context dependence lands.
    """

    delta_i: float = 1e-4       # gap open, insertion
    delta_d: float = 1e-4       # gap open, deletion
    gamma_i: float = 0.3        # gap extend, insertion
    gamma_d: float = 0.3        # gap extend, deletion
    tau: float = 1e-4           # termination

    def __post_init__(self) -> None:
        for f in fields(self):
            v = getattr(self, f.name)
            if not 0.0 <= v < 1.0:
                raise ValueError(f"{f.name} must lie in [0, 1), got {v}")
        if self.delta_i + self.delta_d + self.tau >= 1.0:
            raise ValueError(
                "delta_i + delta_d + tau must be < 1, else M->M is non-positive")
        if self.gamma_i + self.tau >= 1.0:
            raise ValueError(
                "gamma_i + tau must be < 1, else I->M is non-positive")

    def log_trans(self) -> LogTrans:
        return LogTrans(
            b_m=_log(1.0 - self.delta_i),
            b_i=_log(self.delta_i),
            m_m=_log(1.0 - self.delta_i - self.delta_d - self.tau),
            m_i=_log(self.delta_i),
            m_d=_log(self.delta_d),
            m_e=_log(self.tau),
            i_i=_log(self.gamma_i),
            i_m=_log(1.0 - self.gamma_i - self.tau),
            i_e=_log(self.tau),
            d_d=_log(self.gamma_d),
            d_m=_log(1.0 - self.gamma_d),
        )


def eps_from_phred(q):
    """Per-base error probability 10^(-Q/10), clamped to a meaningful range.

    Clamping is not cosmetic: HiFi BAMs use Q0 for "unknown", which would give
    eps = 1 and log(1 - eps) = -inf, poisoning the whole forward pass.
    """
    eps = np.power(10.0, -np.asarray(q, dtype=np.float64) / 10.0)
    return np.clip(eps, EPS_MIN, EPS_MAX)


def mismatch_cost(q):
    """s = log(3(1 - eps)/eps), the cost of one substitution (model.tex eq. 5).

    5.7 nats at Q20, 8.0 at Q30, 10.3 at Q40.  The `exp(-NM)` placeholder in
    `qatool.plotters.likelihood` is s = 1, i.e. eps = 3/(3 + e) = 0.52.
    """
    eps = eps_from_phred(q)
    return np.log(3.0 * (1.0 - eps) / eps)


@dataclass(frozen=True)
class QualityTable:
    """Log emission probabilities indexed by Phred score, 0..MAX_PHRED."""

    log_match: np.ndarray
    log_mismatch: np.ndarray
    log_insert: float

    def cost(self) -> np.ndarray:
        """Per-quality mismatch cost s; equals `mismatch_cost(arange(...))`."""
        return self.log_match - self.log_mismatch


def quality_table(log_insert: float | None = None) -> QualityTable:
    """Precompute emissions for every Phred value.

    `log_insert` defaults to log(1/4), the context-independent choice of §5.4.
    That is crude for HiFi -- inserted bases sit overwhelmingly in homopolymers,
    where the base is near-determined by context -- and is the first thing the
    context-dependent model of §5.7 replaces.
    """
    eps = eps_from_phred(np.arange(MAX_PHRED + 1))
    return QualityTable(
        log_match=np.log1p(-eps),           # accurate for the small eps that dominate
        log_mismatch=np.log(eps / 3.0),
        log_insert=math.log(0.25) if log_insert is None else log_insert,
    )
