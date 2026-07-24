# assemblyQA

Quality-assessment tooling for HiFi genome assemblies. Given a candidate assembly
region and a reference, it recruits the HiFi reads that support each assembly and
produces a diagnostic report plot (k-mer dotplot + read-alignment overlay + marginal
coverage) so the two assemblies can be compared locus by locus.

## Layout

```
src/qatool/          Importable package
  SegmentStorage.py    Segment + SegmentCollection (non-overlapping interval set)
  aln_utils.py         CIGAR parsing / identity + length stats
  build_aln_db.py      BAM -> SQLite `aln` table
  SRRecruiter.py       Recursive read<->segment recruiter + likelihood metrics
  analysis.py          Read-only DB query helpers (return DataFrames)
  plotters.py          Mosaic k-mer dotplot + coverage report plot
  alignment_plot.py    Coverage / alignment plotting helpers
  RecruitLogger.py     Recruitment logging
scripts/             Data-prep shell drivers (prepare_data.sh, run_prepare.sh, align.sh)
docs/
  plotters.md          Deep dev notes on plotters.py (coordinate system, fixed bugs)
  comparative_study.md 52-assembly VGP vs LJA comparative study notes
examples/            Comparative-study code + (git-ignored) data
notebooks (repo root) Interactive harnesses; reference data by paths relative to root
```

Large data/result directories (`alignments/`, `genomes/`, `patching/`, `data/`,
`examples/{genomes,bgz,dbs,...}`, …) live on disk but are excluded from git via
`.gitignore`.

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .          # exposes the `qatool` package
```

`pip install -e .` makes `import qatool` work from anywhere, including the notebooks
at the repo root. Notebooks and scripts import modules as `qatool.<module>`
(e.g. `from qatool.SRRecruiter import Recruiter`).

## Data preparation

```bash
scripts/prepare_data.sh --asm1 <ref.fasta> --asm1-name ref \
                        --asm2 <cand.fasta> --asm2-name cand \
                        --reads <hifi_dir_or_fastq> --label <label> \
                        --outdir <dir> [--threads N]
```

This aligns reads to each assembly (minimap2), bgzip-indexes the FASTAs, and builds a
combined SQLite alignment DB (`qatool.build_aln_db`). The script puts `src/` on
`PYTHONPATH`, so it runs without a prior `pip install`.

## Secrets

The comparative-study sheet updater needs Google OAuth credentials. They are kept in
the untracked file `.secrets/google_oauth.env` (never committed). Load with:

```bash
set -a; source .secrets/google_oauth.env; set +a
```

## Related

`genome_patcher/` (symlink to a separate repo) patches assemblies with curated patch
sequences — see its own README.
