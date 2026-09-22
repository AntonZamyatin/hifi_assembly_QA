"""Pair HMM for the read mark likelihood (model.tex §5).

Score a read against one template position:

    from qamodel.phmm import PhmmParams, quality_table, forward_window
    s = forward_window(read, qual, template, PhmmParams(), quality_table())
    s.logp          # log P(read | this position)

See docs/phmm.md for the staged plan and docs/phmm_explained.pdf for the model.
"""

from qamodel.phmm.band import DEFAULT_WIDTH, Band, cigar_path
from qamodel.phmm.dp import (
    ViterbiResult,
    WindowScore,
    encode,
    forward,
    forward_window,
    viterbi,
)
from qamodel.phmm.kernel import HAVE_NUMBA, forward_numba
from qamodel.phmm.sample import Sample, run_lengths, sample_read, sample_reads
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
    # configuration in
    "PhmmParams", "LogTrans", "QualityTable", "quality_table",
    "Band", "cigar_path", "DEFAULT_WIDTH",
    # scoring
    "forward_numba",    # fastest: compiled, linear space; needs numba
    "forward_window",   # portable: numpy, row-recycled, bandable
    "forward",          # oracle: full table, carries full_template
    "viterbi",          # best single path; never used to score
    "HAVE_NUMBA",
    # results out
    "WindowScore", "ViterbiResult",
    # sampling: the model run forwards
    "sample_read", "sample_reads", "Sample", "run_lengths",
    # helpers
    "encode", "eps_from_phred", "mismatch_cost",
    "EPS_MAX", "EPS_MIN", "MAX_PHRED",
]
