"""Pair HMM for the read mark likelihood (model.tex §5).

See docs/phmm.md for the staged plan.  Stage A: parameters and emissions.
"""

from qamodel.phmm.params import (
    EPS_MAX,
    EPS_MIN,
    MAX_PHRED,
    LogTrans,
    PhmmParams,
    QualityTable,
    eps_from_phred,
    mismatch_cost,
    quality_table,
)

__all__ = [
    "EPS_MAX",
    "EPS_MIN",
    "MAX_PHRED",
    "LogTrans",
    "PhmmParams",
    "QualityTable",
    "eps_from_phred",
    "mismatch_cost",
    "quality_table",
]
