#!/usr/bin/env python3
"""Gather the per-region pair-HMM scores written by the hg002_hic notebooks.

Each region's run_pipe() writes one JSON file to analysis/hg002_hic/phmm_scores/.
This collects them into one TSV with one row per recruited locus, sorted so the
strongest verdicts for either candidate sit at the two ends. Delta S is positive
when the reads favour asm1 (hifiasm_hic) and negative when they favour asm2
(curated).

Seeds inside one repeat family usually recruit the same locus. Counting each
seed as a row would count that locus once per seed, so rows are merged by the
`locus` id run_pipe records (its segments and read set); `regions` lists the
seeds that reached it. JSON from runs before that field existed is grouped by
its segment totals, read counts and Delta S instead, and its id starts with "~".

Usage: venv/bin/python examples/scripts/collect_phmm_scores.py
"""

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

SCORES_DIR = Path("/scratch/projects/assemblyQA/examples/analysis/hg002_hic/phmm_scores")
OUT = SCORES_DIR.parent / "phmm_scores.tsv"
OUT_REGIONS = SCORES_DIR.parent / "phmm_scores_by_region.tsv"   # one row per seed

STATS = ["delta_s", "emission", "exposure", "reads", "dropped",
           "unscorable", "at_plus_C", "at_minus_C", "L1", "L2", "k1", "k2",
           "worst_band_edge", "delta_s_C50", "delta_s_C1000"]
COLUMNS = ["locus", "n_regions", "regions", "title"] + STATS
REGION_COLUMNS = ["region", "locus", "title"] + STATS

def inferred_locus(d: dict) -> str:
    sig = [d["asm1"], d["asm2"], d["L1"], d["L2"], d["k1"], d["k2"], d["reads"],
           d["dropped"], round(d["delta_s"], 4)]
    return "~" + hashlib.sha1(json.dumps(sig).encode()).hexdigest()[:15]


by_locus = defaultdict(list)
for path in sorted(SCORES_DIR.glob("*.json")):
    d = json.loads(path.read_text())
    d["region"] = path.stem
    d["locus"] = d.get("locus") or inferred_locus(d)
    d["delta_s_C50"] = d["sensitivity"]["50"]
    d["delta_s_C1000"] = d["sensitivity"]["1000"]
    by_locus[d["locus"]].append(d)

rows = []
for lid, ds in by_locus.items():
    d = dict(ds[0])
    d["n_regions"] = len(ds)
    d["regions"] = ",".join(x["region"] for x in ds)
    rows.append(d)

def write(path, records, columns):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in records:
            w.writerow({k: (f"{d[k]:.2f}" if isinstance(d[k], float) else d[k]) for k in columns})


rows.sort(key=lambda d: d["delta_s"])
write(OUT, rows, COLUMNS)
write(OUT_REGIONS, sorted((d for ds in by_locus.values() for d in ds), key=lambda d: d["delta_s"]),
      REGION_COLUMNS)

n_regions = sum(d["n_regions"] for d in rows)
fav1 = sum(d["delta_s"] > 0 for d in rows)
fav2 = sum(d["delta_s"] < 0 for d in rows)
flat = len(rows) - fav1 - fav2
print(f"{n_regions} regions in {len(rows)} distinct loci -> {OUT}")
print(f"  one row per region -> {OUT_REGIONS}")
if rows:
    print(f"  favour {rows[0]['asm1']}: {fav1}   favour {rows[0]['asm2']}: {fav2}   exactly 0: {flat}")
    print(f"  C-sensitive (sign differs between C=50 and C=1000): "
          f"{sum((d['delta_s_C50'] > 0) != (d['delta_s_C1000'] > 0) for d in rows)}")
