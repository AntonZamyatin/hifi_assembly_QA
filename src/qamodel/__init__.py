"""qamodel -- implementation of the likelihood model in docs/model.tex.

    from qamodel import score_locus, library_stats, fasta_length
    lib = library_stats(db, "curated", fasta_length("hg002.curated.fasta.bgz"))
    result = score_locus(reads, segments, "hifiasm_hic", "curated",
                         db_path=db, bams=bams, fastas=fastas, library=lib)
    print(result.summary())

Deliberately imports nothing from `qatool`: the dependency runs one way only
(`qatool.plotters.likelihood` -> `qamodel`), which is what keeps the model core
usable on simulated sequence with no database, recruiter or BAM in play.
"""

from qamodel.exposure import Exposure, Library, exposure_term, fasta_length, library_stats
from qamodel.score import DEFAULT_C, LocusScore, ReadScore, score_locus

__all__ = ["score_locus", "LocusScore", "ReadScore", "DEFAULT_C",
           "library_stats", "fasta_length", "Library", "exposure_term", "Exposure"]
