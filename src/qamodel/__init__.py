"""qamodel — implementation of the likelihood model in docs/model.tex.

Deliberately imports nothing from `qatool`: the dependency runs one way only
(`qatool.plotters.likelihood` -> `qamodel`), which is what keeps the model core
usable on simulated sequence with no database, recruiter or BAM in play.
"""
