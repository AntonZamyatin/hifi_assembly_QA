#!/usr/bin/env python3
"""Write qaview project files for the analyses whose notebooks are in examples/analysis.

- hg002: hifiasm against curated, seeds from QUAST (examples/quast/hg002), as
  generate_notebooks_hg002.py.
- one project per species of generate_notebooks.py's IG batch: ref against lja,
  seeds from anomaly_table.csv (every IG locus of that species), written as a
  BED file each under projects/seeds/ and named like the notebook figures,
  `<locus>_<i:03d>_<contig>_<type>` with `i` counting rows of one (species, locus).
  A row that is a single position (min_coord = max_coord; the notebooks passed
  it as a zero-length segment) becomes the 1 bp interval at that position.

- figure_for_anton: mPleAur1 reference against candidate 2, the seed of
  "Figure for Anton.ipynb" in the repository root (its data are under the root's
  alignments/ and genomes/, not examples/).

Recruitment is the hg002_hic project's [recruitment] block with the anchor
extension; a species whose BAMs or FASTAs are missing is skipped with a note.

Run: venv/bin/python examples/scripts/make_projects.py
"""

import csv
import os
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]          # examples/
PROJECTS = BASE / "projects"
SEEDS = PROJECTS / "seeds"

# The IG batch of generate_notebooks.py (OUR_LABELS) and its anomaly table.
IG_LABELS = ["mApoSyl1", "mCanLor1", "mCynVol1", "mDasNov1", "mDelDel1", "mEptNil1",
             "mMyoDau2", "mMyoYum1", "mOryCun1", "mPleAur1", "rEryReg1", "rVipUrs1"]
ANOMALY_CSV = BASE / "anomaly_table.csv"

RECRUITMENT = """\
[recruitment]                             # hg002_hic's values, with the anchor extension
pid_threshold    = 0.997
plen_threshold   = 0.999
best_pid_delta   = 0.001
merge_gap        = 50000
extension_window = 100
ext_max          = 300000
ext_rounds       = 100
max_rounds       = 100
extension        = "anchor"

[resources]
recruit_workers = 32
score_workers   = 32
"""


def candidate(key: str, name: str, label: str, note: str = "") -> str:
    head = f"[{key}]"
    return (f"{head:<42}{'# ' + note if note else ''}".rstrip() + "\n"
            f'name  = "{name}"\n'
            f'fasta = "../bgz/{label}.{name}.fasta.bgz"\n'
            f'bam   = "../alignments/{label}.{name}.sorted.bam"\n')


def missing_files(label: str, names) -> list:
    want = [BASE / "dbs" / f"{label}.aln.sqlite"]
    for n in names:
        want += [BASE / "bgz" / f"{label}.{n}.fasta.bgz", BASE / "bgz" / f"{label}.{n}.fasta.bgz.fai",
                 BASE / "alignments" / f"{label}.{n}.sorted.bam",
                 BASE / "alignments" / f"{label}.{n}.sorted.bam.bai"]
    return [str(p.relative_to(BASE)) for p in want if not p.exists()]


def write(path: Path, text: str):
    path.write_text(text)
    print(f"wrote {path.relative_to(BASE)}")


def hg002():
    miss = missing_files("hg002", ("hifiasm", "curated"))
    if miss:
        print(f"skipped hg002: missing {', '.join(miss)}")
        return
    write(PROJECTS / "hg002.toml", f"""\
# qaview project: HG002 hifiasm against the curated assembly, the analysis of
# examples/analysis/hg002 (generate_notebooks_hg002.py). Written by
# examples/scripts/make_projects.py. Paths are relative to this file.

name = "hg002"

{candidate("asm1", "hifiasm", "hg002", "Delta S > 0 favours asm1")}
{candidate("asm2", "curated", "hg002")}
[data]
db = "../dbs/hg002.aln.sqlite"
library_assembly = "curated"              # lambda reference; keep fixed across loci

[[seeds]]
format   = "quast"
path     = "../quast/hg002/misassemblies.tsv"
assembly = "hifiasm"

{RECRUITMENT}
""")


def ig():
    rows = defaultdict(list)
    with open(ANOMALY_CSV) as fh:
        for r in csv.DictReader(fh):
            if r["assembly"] in IG_LABELS:
                rows[r["assembly"]].append(r)
    SEEDS.mkdir(parents=True, exist_ok=True)
    for label in IG_LABELS:
        miss = missing_files(label, ("ref", "lja"))
        if miss:
            print(f"skipped {label}: missing {', '.join(miss)}")
            continue
        count = defaultdict(int)
        lines = [f"# {label}: IG anomalies from anomaly_table.csv, on ref (make_projects.py)"]
        for r in rows[label]:
            count[r["locus"]] += 1
            name = f"{r['locus']}_{count[r['locus']]:03d}_{r['contig']}_{r['type']}"
            start, end = int(r["min_coord"]), int(r["max_coord"])
            end = max(end, start + 1)               # a single position: its 1 bp
            lines.append(f"{r['contig']}\t{start}\t{end}\t{name}")
        bed = SEEDS / f"{label}_ig.bed"
        write(bed, "\n".join(lines) + "\n")
        loci = ", ".join(sorted(count))
        write(PROJECTS / f"{label}_ig.toml", f"""\
# qaview project: {label} reference against its LJA assembly over the IG loci
# ({loci}), the analysis of examples/analysis/{label}_*.ipynb (generate_notebooks.py).
# Written by examples/scripts/make_projects.py. Paths are relative to this file.

name = "{label}_ig"

{candidate("asm1", "ref", label, "Delta S > 0 favours asm1")}
{candidate("asm2", "lja", label)}
[data]
db = "../dbs/{label}.aln.sqlite"
library_assembly = "ref"                  # lambda reference; keep fixed across loci

[[seeds]]
format   = "bed"
path     = "seeds/{label}_ig.bed"
assembly = "ref"

{RECRUITMENT}
""")


def figure_for_anton():
    root = BASE.parent
    files = {"db": root / "alignments" / "mPleAur1.aln.sqlite",
             "ref_bam": root / "alignments" / "mPleAur1.ref.sorted.bam",
             "cand2_bam": root / "alignments" / "mPleAur1.cand2.sorted.bam",
             "ref_fasta": root / "genomes" / "mPleAur1" / "mPleAur1.dip.ref.fasta.bgz",
             "cand2_fasta": root / "genomes" / "mPleAur1" / "mPleAur1.dip.cand2.fasta.bgz"}
    miss = [str(f.relative_to(root)) for f in files.values() if not f.exists()]
    if miss:
        print(f"skipped figure_for_anton: missing {', '.join(miss)}")
        return
    SEEDS.mkdir(parents=True, exist_ok=True)
    write(SEEDS / "figure_for_anton.bed",
          "# the seed of 'Figure for Anton.ipynb' (repository root), on ref\n"
          "SUPER_1\t1800000\t1900000\tSUPER_1_1800000_1900000\n")
    rel = lambda f: os.path.relpath(f, PROJECTS)
    write(PROJECTS / "figure_for_anton.toml", f"""\
# qaview project: mPleAur1 reference against candidate 2 at the seed of
# "Figure for Anton.ipynb" in the repository root, whose data are under the root's
# alignments/ and genomes/. Written by examples/scripts/make_projects.py.
# Paths are relative to this file.

name = "figure_for_anton"

[asm1]                                    # Delta S > 0 favours asm1
name  = "ref"
fasta = "{rel(files['ref_fasta'])}"
bam   = "{rel(files['ref_bam'])}"

[asm2]
name  = "cand2"
fasta = "{rel(files['cand2_fasta'])}"
bam   = "{rel(files['cand2_bam'])}"

[data]
db = "{rel(files['db'])}"
library_assembly = "ref"                  # lambda reference; keep fixed across loci

[[seeds]]
format   = "bed"
path     = "seeds/figure_for_anton.bed"
assembly = "ref"

{RECRUITMENT}
""")


if __name__ == "__main__":
    hg002()
    ig()
    figure_for_anton()
