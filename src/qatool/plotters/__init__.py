"""Plotting for the HiFi assembly QA report.

Public API is re-exported here so existing imports keep working:
    from qatool.plotters import report_plot, alignment_plot, ...
The implementation is split across submodules (kmer, coverage, reads_overlay,
classify, likelihood, report, axis, db, colors).
"""

from qatool.plotters.colors import (
    SUPPORT_COLORS,
    OCCURRENCE_COLORS,
    PID_TIE_EPS_DEFAULT,
    BEST_PID_DELTA_DEFAULT,
)
from qatool.plotters.kmer import dotplot_kmer, dotplot_kmer_mosaic
from qatool.plotters.coverage import alignment_plot
from qatool.plotters.reads_overlay import reads_alignment_dotplot
from qatool.plotters.likelihood import likelihood_report
from qatool.plotters.report import report_plot

__all__ = [
    "report_plot",
    "alignment_plot",
    "reads_alignment_dotplot",
    "dotplot_kmer",
    "dotplot_kmer_mosaic",
    "likelihood_report",
    "SUPPORT_COLORS",
    "OCCURRENCE_COLORS",
    "PID_TIE_EPS_DEFAULT",
    "BEST_PID_DELTA_DEFAULT",
]
