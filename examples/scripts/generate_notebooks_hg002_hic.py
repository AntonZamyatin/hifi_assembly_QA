#!/usr/bin/env python3
"""Generate one Jupyter notebook per (contig, misassembly_type) for hg002 hifiasm_hic assembly."""

import csv
import json
import uuid
from collections import defaultdict
from pathlib import Path

BASE         = Path("/scratch/projects/assemblyQA/examples")
CODE_DIR     = Path("/scratch/projects/assemblyQA")
ANALYSIS_DIR = BASE / "analysis" / "hg002_hic"
PICS_DIR     = ANALYSIS_DIR / "pics"
SCORES_DIR   = ANALYSIS_DIR / "phmm_scores"
LOCI_DIR     = ANALYSIS_DIR / "loci"      # pair-HMM scores keyed by locus content
MISASM_TSV   = BASE / "quast" / "hg002_hic" / "misassemblies.tsv"
PIC_PREFIX   = "reads_"   # distinguishes new reads-lines figures from the old k-mer ones

IMPORTS = """\
import json
import os

from IPython.display import Image, display
import matplotlib.pyplot as plt
import pysam
import pandas as pd
import numpy as np

from qatool.SegmentStorage import *
from qatool.SRRecruiter import Recruiter, metrics_report_from_aln_df
from qatool.analysis import *
from qatool.alignment_plot import *
from qatool.plotters import *
"""

ANALYSIS_PARAMS = """\
# ── analysis parameters (edit here to tune the pipeline) ─────────────────────
pid_threshold      = 0.997   # recruitment
plot_pid_threshold = 0.997   # plotting / classification (match recruitment)
plen_threshold     = 0.999
best_pid_delta = 0.001
merge_gap      = 50000
ext_window     = 100
ext_max        = 300000
ext_rounds     = 100
max_rounds     = 100
pid_tie_eps    = 0.001
# binsize: automatic (report_plot computes it from segment lengths)
"""

RUN_PIPE = """\
def run_pipe(db_path, asm1, asm2, bam1, bam2, fasta1, fasta2, sg, fig_path,
             title=None, score_path=None):
    SRR = Recruiter(
        db_path=db_path, asm1=asm1, asm2=asm2,
        asm1_bam_path=bam1,
        asm2_bam_path=bam2,
        asm1_fasta_path=fasta1,
        asm2_fasta_path=fasta2,
        merge_gap=merge_gap,
        verbose=True,
        pid_threshold=pid_threshold, plen_threshold=plen_threshold,
        best_pid_delta=best_pid_delta,
        extension_window=ext_window, ext_max=ext_max, ext_rounds=ext_rounds,
        max_rounds=max_rounds, out_dir=None,   # per-round plots are disabled
    )
    segments, reads, df = SRR.run(sg)

    ref_seq, cand_seq = get_concatenated_segments(
        sorted(segments, key=lambda s: (s.contig, s.start, s.end)),
        fasta1, fasta2,
        asm1=asm1, asm2=asm2,
    )

    report_plot(
        ref_seq, cand_seq, segments, reads,
        db_path, plot_pid_threshold, plen_threshold,
        asm1_bam=bam1,
        asm2_bam=bam2,
        pid_tie_eps=pid_tie_eps,
        best_pid_delta=best_pid_delta,
        boundary_params=("black", 1.0, 1),
        asm1=asm1, asm2=asm2, asm1_fasta=fasta1, asm2_fasta=fasta2,
        support_stripes=True,
        classification="support",
        output_path=fig_path,
        title=title,
        show=False,
    )

    # pair-HMM locus score: emission (clamped at +/-C) + exposure; positive favours asm1.
    # Seeds inside one repeat family usually recruit the same locus; the score is
    # cached by locus content in LOCI_DIR, so each distinct locus is scored once.
    result = phmm_likelihood_report(
        segments, reads, asm1, asm2,
        db_path=db_path, bams={asm1: bam1, asm2: bam2},
        fastas={asm1: fasta1, asm2: fasta2},
        genome_fasta=genome_fasta, library_assembly=library_assembly,
        cache_dir=LOCI_DIR,
    )
    if score_path:
        up, down = result.clamped_counts()
        x = result.exposure
        with open(score_path, "w") as fh:
            json.dump({
                "title": title, "asm1": asm1, "asm2": asm2,
                "locus": locus_id(segments, reads, asm1, asm2),
                "delta_s": result.delta_s, "emission": result.emission,
                "exposure": x.value, "C": result.C,
                "reads": len(result.reads), "dropped": len(result.dropped),
                "unscorable": len(result.unscorable),
                "at_plus_C": up, "at_minus_C": down,
                "worst_band_edge": result.worst_band_edge,
                "L1": x.L1, "L2": x.L2, "k1": x.k1, "k2": x.k2,
                "sensitivity": {str(c): result.at(c) for c in (50, 100, 200, 500, 1000)},
            }, fh, indent=1)
    return result
"""


def make_id():
    return str(uuid.uuid4()).replace("-", "")[:16]


def code_cell(source):
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": make_id(),
        "metadata": {},
        "outputs": [],
        "source": source,
    }


def md_cell(source):
    return {
        "cell_type": "markdown",
        "id": make_id(),
        "metadata": {},
        "source": source,
    }


def make_notebook(cells):
    return {
        "cells": cells,
        "metadata": {
            # the project venv, registered with
            #   venv/bin/python -m ipykernel install --user --name assemblyqa
            # the stock python3 kernel is mambaforge base and cannot import qatool
            "kernelspec": {
                "display_name": "Python (assemblyQA venv)",
                "language": "python",
                "name": "assemblyqa",
            },
            "language_info": {
                "name": "python",
                "version": "3.10.0",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


# ── load and group misassemblies ──────────────────────────────────────────────
groups = defaultdict(list)
with open(MISASM_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        groups[(row["contig"], row["type"])].append(row)

ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
PICS_DIR.mkdir(parents=True, exist_ok=True)
SCORES_DIR.mkdir(parents=True, exist_ok=True)
LOCI_DIR.mkdir(parents=True, exist_ok=True)

for (contig, mis_type), rows in sorted(groups.items()):
    nb_stem = f"{contig}_{mis_type}"

    common_params = f"""\
BASE  = "{BASE}"
LABEL = "hg002"
PICS_DIR = "{PICS_DIR}"

db_path = f"{{BASE}}/dbs/{{LABEL}}.aln.sqlite"
asm1    = "hifiasm_hic"
asm2    = "curated"
bam1    = f"{{BASE}}/alignments/{{LABEL}}.hifiasm_hic.sorted.bam"
bam2    = f"{{BASE}}/alignments/{{LABEL}}.curated.sorted.bam"
fasta1  = f"{{BASE}}/bgz/{{LABEL}}.hifiasm_hic.fasta.bgz"
fasta2  = f"{{BASE}}/bgz/{{LABEL}}.curated.fasta.bgz"
SCORES_DIR = "{SCORES_DIR}"
LOCI_DIR   = "{LOCI_DIR}"

# lambda for the pair-HMM score: primary reads on the curated diploid assembly.
# Keep these fixed across every locus -- lambda must not depend on the locus.
genome_fasta     = fasta2
library_assembly = "curated"
"""

    header = (
        f"# hg002 hifiasm_hic — {mis_type} on `{contig}`\n\n"
        f"**Misassemblies:** {len(rows)}"
    )

    cells = [
        md_cell(header),
        code_cell(IMPORTS),
        code_cell(ANALYSIS_PARAMS),
        code_cell(RUN_PIPE),
        code_cell(common_params),
    ]

    for i, r in enumerate(rows, 1):
        fig_name    = f"{PIC_PREFIX}{nb_stem}_{i:03d}.png"
        bp_start = int(r['bp_start'])
        bp_end   = int(r['bp_end'])
        bp_str   = f"{bp_start:,}" if bp_start == bp_end else f"{bp_start:,} – {bp_end:,}"
        incon = f"  \n**Inconsistency:** {r['inconsistency_bp']} bp" if r["inconsistency_bp"] else ""
        region_header = f"## Region {i}: breakpoint {bp_str} (window {int(r['start']):,} – {int(r['end']):,}){incon}"
        plot_title = f"{r['contig']} — {mis_type} — bp {bp_str}"
        region_cell = f"""\
sg = SegmentCollection([Segment("hifiasm_hic", "{r['contig']}", {r['start']}, {r['end']})])
run_pipe(db_path, asm1, asm2, bam1, bam2, fasta1, fasta2, sg, f"{{PICS_DIR}}/{fig_name}",
         title="{plot_title}", score_path=f"{{SCORES_DIR}}/{nb_stem}_{i:03d}.json")"""

        cells.append(md_cell(region_header))
        cells.append(code_cell(region_cell))

    nb_path = ANALYSIS_DIR / f"{nb_stem}.ipynb"
    with open(nb_path, "w") as f:
        json.dump(make_notebook(cells), f, indent=1)
    print(f"  {nb_path.name}  ({len(rows)} regions)")

print(f"\nGenerated {len(groups)} notebooks → {ANALYSIS_DIR}")
