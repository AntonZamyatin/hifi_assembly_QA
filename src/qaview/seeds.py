"""Seed regions from a QUAST misassembly table, for picking in the sidebar.

The table is the one `examples/quast/extract_misassemblies.py` writes and the
analysis notebooks are generated from: one row per misassembly on an asm1 contig,
with the +/-20 kb window used as the recruitment seed and the breakpoint inside it.

Each row gets the same region name the notebooks use, `<contig>_<type>_<i:03d>`
with `i` counting rows of one (contig, type) in file order, so a seed here can be
matched to its notebook figure and to the pair-HMM score the notebook wrote.

Conventional locations, relative to the data `base` and a QUAST run name:

    <base>/quast/<run>/misassemblies.tsv
    <base>/analysis/<run>/phmm_scores/<region>.json     (optional)
"""

import csv
import json
import os
from collections import defaultdict

import pandas as pd

COLUMNS = ["seed", "region", "contig", "type", "breakpoint", "bp", "start", "end",
           "inconsistency_bp", "delta_s"]


# Short type names for the sidebar's narrow `seed` column.
ABBREV = {"translocation": "transl", "relocation": "reloc", "inversion": "inv",
          "interspecies translocation": "inter"}


def seeds_path(base: str, run: str) -> str:
    return os.path.join(base, "quast", run, "misassemblies.tsv")


def scores_dir(base: str, run: str) -> str:
    return os.path.join(base, "analysis", run, "phmm_scores")


def _delta_s(scores: str | None, region: str):
    if not scores:
        return None
    path = os.path.join(scores, f"{region}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return round(json.load(fh)["delta_s"], 2)


def load_seeds(path: str, scores: str | None = None) -> pd.DataFrame:
    """The misassembly table as seeds, one row each, in file order.

    `delta_s` is the notebook's pair-HMM score for the region when `scores` holds
    one, else missing.  Those scores used the notebooks' recruitment parameters,
    which need not match the sidebar's.
    """
    rows, seen = [], defaultdict(int)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["contig"], r["type"])
            seen[key] += 1
            region = f"{r['contig']}_{r['type']}_{seen[key]:03d}"
            bp0, bp1 = int(r["bp_start"]), int(r["bp_end"])
            rows.append(dict(
                seed=f"{r['contig']} {ABBREV.get(r['type'], r['type'])}{seen[key]}",
                region=region, contig=r["contig"], type=r["type"],
                breakpoint=f"{bp0:,}" if bp1 - bp0 <= 1 else f"{bp0:,}-{bp1:,}",
                bp=bp0,
                start=int(r["start"]), end=int(r["end"]),
                inconsistency_bp=int(r["inconsistency_bp"]) if r["inconsistency_bp"] else None,
                delta_s=_delta_s(scores, region),
            ))
    return pd.DataFrame(rows, columns=COLUMNS)
