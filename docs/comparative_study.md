# Assembly QA Project — Session Context

**Last updated:** 2026-06-03  
**Primary working directory:** `/scratch/projects/assemblyQA/examples/`

---

## Project Overview

Comparative genomics assembly quality assessment study. The goal is to assess assembly anomalies (immunoglobulin locus breaks, mismatches) across a set of 52 VGP mammal + reptile assemblies by comparing VGP curated assemblies against LJA de-novo assemblies assembled from the same HiFi reads.

**Google Sheet (master tracking table):**  
`https://docs.google.com/spreadsheets/d/1_Sok7nzTtYpAsIH9BYLlHvjMyWe5eNMk5kHx8O9caj4`  
Sheet2 GID: `447612736`

---

## Directory Layout

```
/scratch/projects/assemblyQA/examples/
│                                 (documented in ../docs/comparative_study.md)
├── anomaly_table.csv             ← 117-row table: assembly,contig,min_coord,max_coord,locus,type
├── anomaly_table_with_haplome.tsv
├── find_anomaly_haplomes.py
├── setup_genomes.py
├── scripts/
│   ├── update_refs_and_contigs.py   ← updates Sheet2 cols C/D (ref paths) and I (contigs found)
│   ├── update_sheet.py              ← earlier sheet updater
│   ├── download_cur_assemblies.py   ← downloads VGP curated assemblies
│   ├── download_hifi.py             ← downloads HiFi reads from SRA/ENA
│   ├── download_vgp_assemblies.py   ← downloads VGP assemblies
│   └── query_hifi_sra.py            ← queries SRA for HiFi run metadata
├── tables/
│   ├── hifi_sra_table.tsv           ← 48 rows: vgp_id, sci_name, hifi_data_dir, biosample, sra_runs, ftp_urls, total_gb, status
│   └── hifi_status.csv              ← 52 rows: per-assembly HiFi read status + assembly status
├── genomes/
│   ├── {vgp_id}.{hap}.fasta.gz      ← symlinks to curated VGP assemblies (86 files)
│   └── {vgp_id}.lja.fasta           ← symlinks to LJA assemblies (14 files as of 2026-06-03)
└── lja/
    └── {label}/assembly.fasta       ← LJA assembly output (flat, NOT experimental/ subdir)
```

**Mass assembly script:**  
`/scratch/projects/assembly_jobs/run_lja_all.sh`  
33 jobs, 5 parallel, 40 threads/job. Skip condition: `lja/{label}/assembly.fasta` exists.

---

## Naming Conventions

### VGP genome symlinks in `genomes/`
Pattern: `{vgp_id}.{hap_tag}.fasta.gz`

Hap tags:
- **ref1** (haplotype 1): `pri`, `hap1`, `mat`
- **ref2** (haplotype 2): `alt`, `hap2`, `pat`

Original filenames follow: `{vgp_id}.{hap_tag}.cur.{YYYYMMDD}.fasta.gz`  
Some NCBI assemblies are `.fna.gz` — these use NCBI accession IDs in headers with VGP names in descriptions. The script mines descriptions with regex `(?<!\w)((?:SUPER|SCAFFOLD|scaffold|atg|hap_ptg|ptg|ctg|h\dtg)[^\s,;]+)`.

### LJA symlinks in `genomes/`
Pattern: `{vgp_id}.lja.fasta` → `lja/{label}/assembly.fasta`  
Exception: `mMyoYum1.lja.fasta` → `lja/mMyoYum1.0/assembly.fasta` (directory has `.0` suffix)

---

## 52-Assembly Tracking List (Sheet2 row order)

| Row | VGP ID    | HiFi status              | LJA status (2026-06-03)     |
|-----|-----------|--------------------------|------------------------------|
| 2   | mApoSyl1  | raw files (4 SRR, 145G)  | **DONE**                     |
| 3   | mCanLor1  | raw files (7 SRR, 366G)  | **DONE**                     |
| 4   | mCanLor2  | no data                  | —                            |
| 5   | mChiNiv1  | merged (43G)             | pending                      |
| 6   | mCorTow1  | merged (52G)             | pending                      |
| 7   | mCynVol1  | raw uncompressed (210G)  | **DONE**                     |
| 8   | mDasNov1  | raw files (5 SRR, 75G)   | **DONE**                     |
| 9   | mDelDel1  | raw files (4 SRR, 190G)  | **DONE**                     |
| 10  | mDicBic1  | merged (618G)            | pending                      |
| 11  | mEptNil1  | raw files (4 SRR, 125G)  | **DONE**                     |
| 12  | mEriEur2  | merged (60G)             | pending                      |
| 13  | mGloMel1  | merged (44G)             | pending                      |
| 14  | mHetBru1  | merged (83G)             | pending                      |
| 15  | mHipAmp2  | need download (~250G)    | pending                      |
| 16  | mHypAmp2  | merged (51G)             | pending                      |
| 17  | mLagAlb1  | merged (54G)             | pending                      |
| 18  | mLemCat1  | merged (497G)            | pending                      |
| 19  | mMacEug1  | merged (479G)            | pending                      |
| 20  | mManPen7  | merged (269G)            | pending                      |
| 21  | mMarMar1  | merged (43G)             | pending                      |
| 22  | mMicCal1  | not in SRA               | —                            |
| 23  | mMicMin1  | merged (40G)             | pending                      |
| 24  | mMirAng1  | merged (63G)             | pending                      |
| 25  | mMonDom1  | merged (479G)            | pending                      |
| 26  | mMusAve1  | merged (59G)             | pending                      |
| 27  | mMusNiv1  | need download (~268G)    | —                            |
| 28  | mMyoDau2  | raw files (4 SRR, 160G)  | **DONE**                     |
| 29  | mMyoYum1  | raw files (1 SRR, 143G)  | **DONE** (dir: mMyoYum1.0)   |
| 30  | mNeoNeb1  | not in SRA               | —                            |
| 31  | mNycCou1  | not in SRA               | —                            |
| 32  | mOryCun1  | raw files (2 SRR, 173G)  | **DONE**                     |
| 33  | mPerMan1  | merged (60G)             | pending                      |
| 34  | mPhoPho1  | merged (29G)             | pending                      |
| 35  | mPleAur1  | raw files (6 SRR, 526G)  | **DONE**                     |
| 36  | mSorAra2  | merged (194G, subdir)    | pending                      |
| 37  | mSteCoe1  | merged (56G)             | pending                      |
| 38  | mTalEur1  | merged (31G)             | pending                      |
| 39  | mUrsAme1  | merged (57G)             | pending                      |
| 40  | mUrsArc2  | no data                  | —                            |
| 41  | mVesMur1  | merged (40G)             | pending                      |
| 42  | rAllMis2  | no data                  | —                            |
| 43  | rCarCar2  | merged (409G)            | pending                      |
| 44  | rEmyOrb1  | merged (294G)            | **DONE**                     |
| 45  | rEryReg1  | need download (~219G)    | —                            |
| 46  | rLiaOli1  | merged (187G)            | pending                      |
| 47  | rLiaOli2  | no data                  | —                            |
| 48  | rMalTer1  | merged (241G)            | pending                      |
| 49  | rPodCre2  | merged (30G)             | pending                      |
| 50  | rRhiFlo1  | need download (~186G)    | pending (raw SRR files)      |
| 51  | rVipLat1  | merged (559G)            | pending                      |
| 52  | rVipUrs1  | merged (31G)             | **DONE**                     |
| 53  | rZooViv1  | merged (58G)             | pending                      |

**mSorAra2 reads location:** `/data/mammals/sorex_araneus/mSorAra2/hifi_data/` (note subdirectory)  
**mSorAra1** (different individual, already assembled) is in `lja/mSorAra1/` — NOT in the 52-list.

---

## HiFi Data Locations

Base: `/data/mammals/{species_name}/hifi_data/` for mammals  
Base: `/data/reptiles/{species_name}/hifi_data/` for reptiles  
Exceptions:
- `mCanLor1` → `/data/mammals/mCanLor1/hifi_data/`
- `mSorAra2` → `/data/mammals/sorex_araneus/mSorAra2/hifi_data/`

Read selection priority (implemented in `run_lja_all.sh`):
1. `{label}.hifi.fastq.gz` (merged, gzipped)
2. `{label}.hifi.fastq` (merged, uncompressed — mCynVol1, mCanLor1)
3. All `*.fastq` / `*.fastq.gz` in dir (raw SRR runs)
- `.fasta` files excluded (processed artefacts like `.hpc.fasta`)

---

## Google Sheets API Credentials

OAuth2 credentials are **not** stored in version control. They live in the untracked
file `.secrets/google_oauth.env` at the repo root (`SPREADSHEET_ID`, `SHEET_GID`,
`CLIENT_ID`, `REFRESH_TOKEN`). Load them before running the sheet-updater scripts:

```bash
set -a; source .secrets/google_oauth.env; set +a
```

> The scripts under `examples/scripts/` (e.g. `update_refs_and_contigs.py`) still
> contain these values hard-coded — read them from the environment instead and delete
> the literals. The exposed refresh token should be rotated.

Sheet2 column layout:
- A: vgp_id
- B: species name
- C: ref1 path (hap1/pri/mat)
- D: ref2 path (hap2/alt/pat)
- E: lja path (not yet filled for new completions)
- I: contigs found (X/Y format, from anomaly_table.csv scan)

---

## Anomaly Table

`anomaly_table.csv` — 117 rows covering assemblies with detected IG-locus anomalies.  
Columns: `assembly,contig,min_coord,max_coord,locus,type`  
Types: `break`, `mismatch`

Known contig name issues resolved in previous sessions:
- NCBI `.fna.gz` files use accession IDs in headers; VGP names are in description field → mined via DESC_TOKEN regex
- Some assemblies need older pre-curation versions to match anomaly table contig names (mMusNiv1, mOryCun1 use older scaffold naming; mPerMan1 uses GCA v1 accessions)
- Lookup fallbacks: as-is → strip `alt_` prefix → strip leading `_`

---

## Pending Tasks

1. **Run mass LJA assembly:**  
   `nohup bash /scratch/projects/assembly_jobs/run_lja_all.sh &`  
   After completion: create `genomes/{label}.lja.fasta` symlinks for new assemblies and update Sheet2 col E.

2. **Download missing HiFi reads:**
   - `mHipAmp2` — ~250G (download script ready)
   - `mMusNiv1` — ~268G
   - `rEryReg1` — ~219G
   - `rRhiFlo1` — 9 raw SRR files present but needs merging check

3. **Species with no SRA data** (not in 52-list scope yet):
   - `mMicCal1`, `mNeoNeb1`, `mNycCou1` — not in SRA

4. **No data yet** (deferred):
   - `mCanLor2`, `mUrsArc2`, `rAllMis2`, `rLiaOli2`

5. **Update Sheet2 col E (lja paths)** for completed assemblies once runs finish.
