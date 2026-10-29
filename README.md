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
  SRRecruiter.py       Recursive read<->segment recruiter
  plotters/            Mosaic k-mer dotplot + coverage report plot
src/qaview/          Interactive viewer (`python -m qaview project.toml`)
src/qamodel/         Pair-HMM read scoring
scripts/             Data-prep shell drivers (prepare_data.sh, run_prepare.sh)
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

One conda environment holds the Python packages and the external tools that the data
preparation uses (minimap2, samtools, bgzip):

```bash
git clone git@github.com:AntonZamyatin/hifi_assembly_QA.git && cd hifi_assembly_QA
git checkout v0.1.0
mamba env create -f environment.yml      # or conda; once, a few minutes
mamba activate assemblyqa
```

Without conda, `pip install .` (Python >= 3.11) installs everything but the external
tools; the viewer does not need them.

## Run qaview

```bash
qaview <project.toml> [--port 5006] [--cache <dir>]
```

With no project given, the page asks for one.

**From a laptop.** `scripts/qaview-connect.sh` starts qaview on the server and tunnels
to it in one command; the laptop needs only bash, ssh and curl (macOS, Linux, WSL or
Git Bash), not this repository:

```bash
curl -O https://raw.githubusercontent.com/AntonZamyatin/hifi_assembly_QA/master/scripts/qaview-connect.sh
bash qaview-connect.sh <host> --port 5011 --open
```

The page opens on a form to choose a project (`--project <file>` opens one directly).

`<host>` is anything ssh accepts. It finds the clone's `venv` or the conda env
`assemblyqa` on the server (`--python` otherwise); users sharing a server each pick a
`--port`. Ctrl-C stops qaview there. Defaults can live in `QAVIEW_HOST`, `QAVIEW_PORT`,
`QAVIEW_DIR` and `QAVIEW_PYTHON`. By hand instead: `ssh -L 5006:localhost:5006 <host>`,
run `qaview` there, open http://localhost:5006.

**Results go to the cache.** Recruitment, scores, seeds added in the page and
remembered views of a project are kept in `<cache_dir>/<project name>`. `cache_dir` is
set in `qaview_config.toml` beside qaview (default `cache/` in this clone; `--cache`
overrides it for one run). Each clone has its own cache, so another user's project
opens where it lies -- its file names the data by paths relative to itself -- and
only your cache is written:

```bash
qaview /scratch/projects/assemblyQA/examples/projects/ferrets.toml
```

To start from the owner's finished results rather than recomputing, copy them first:
`cp -r /scratch/projects/assemblyQA/cache/ferrets cache/`.

## Development install

```bash
python3 -m venv venv && venv/bin/pip install -e .
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
`PYTHONPATH`, so it runs without a prior `pip install`. The tools come from `PATH` (the
`assemblyqa` environment); set `MINIMAP2`, `SAMTOOLS`, `BGZIP` or `PYTHON` to override.
