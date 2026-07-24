#!/bin/bash
set -euo pipefail

# ── tool paths ───────────────────────────────────────────────────────────────
MINIMAP2=/scratch/azz5280/mambaforge/bin/minimap2
SAMTOOLS=samtools
BGZIP=bgzip
PYTHON=/scratch/azz5280/mambaforge/bin/python3

THREADS=40
OUTDIR=.
ALN_DIR=""
DB_DIR=""
BGZ_DIR=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --asm1 <fasta> --asm1-name <name> \\
                         --asm2 <fasta> --asm2-name <name> \\
                         --reads <fastq_or_dir> --label <label> \\
                         [--outdir <dir>] [--aln-dir <dir>] [--db-dir <dir>] [--bgz-dir <dir>] \\
                         [--threads <n>]

Outputs:
  <aln-dir>/<label>.<name>.sorted.bam[.bai]     sorted HiFi alignments
  <bgz-dir>/<label>.<name>.fasta.bgz[.fai/.gzi] bgzip-indexed FASTA
  <db-dir>/<label>.aln.sqlite                   combined alignment database

--aln-dir / --db-dir / --bgz-dir default to --outdir if not set.
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --asm1)      ASM1_FASTA="$2"; shift 2 ;;
        --asm1-name) ASM1_NAME="$2";  shift 2 ;;
        --asm2)      ASM2_FASTA="$2"; shift 2 ;;
        --asm2-name) ASM2_NAME="$2";  shift 2 ;;
        --reads)     READS="$2";      shift 2 ;;
        --label)     LABEL="$2";      shift 2 ;;
        --outdir)    OUTDIR="$2";     shift 2 ;;
        --aln-dir)   ALN_DIR="$2";    shift 2 ;;
        --db-dir)    DB_DIR="$2";     shift 2 ;;
        --bgz-dir)   BGZ_DIR="$2";    shift 2 ;;
        --threads)   THREADS="$2";    shift 2 ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

[[ -z "${ASM1_FASTA:-}" || -z "${ASM1_NAME:-}" || \
   -z "${ASM2_FASTA:-}" || -z "${ASM2_NAME:-}" || \
   -z "${READS:-}"      || -z "${LABEL:-}" ]] && usage

ALN_DIR="${ALN_DIR:-$OUTDIR}"
DB_DIR="${DB_DIR:-$OUTDIR}"
BGZ_DIR="${BGZ_DIR:-$OUTDIR}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
export PYTHONPATH="$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$ALN_DIR" "$DB_DIR" "$BGZ_DIR"

# ── collect reads (file or directory) ────────────────────────────────────────
collect_reads() {
    local reads=()
    if [[ -f "$READS" ]]; then
        reads+=("$READS")
    else
        while IFS= read -r -d '' f; do
            reads+=("$f")
        done < <(find "$READS" -maxdepth 1 -type f \( \
            -name "*.fasta" -o -name "*.fasta.gz" -o \
            -name "*.fastq" -o -name "*.fastq.gz" \) -print0 | sort -z)
    fi
    [[ ${#reads[@]} -eq 0 ]] && { echo "ERROR: no reads (.fasta/.fastq/.gz) found in $READS"; exit 1; }
    echo "${reads[@]}"
}

# ── step 1: align reads to each assembly (skip if BAM + index exist) ─────────
align_one() {
    local fasta="$1"
    local key="$2"    # e.g. mPleAur1.ref
    local bam="$ALN_DIR/${key}.sorted.bam"

    if [[ -f "$bam" && -f "${bam}.bai" ]]; then
        echo "[align] skip ${key}: $bam already exists"
        return
    fi

    read -ra reads <<< "$(collect_reads)"

    echo "[align] ${key}: $(echo "${reads[@]}" | xargs -n1 basename | tr '\n' ' ')"
    $MINIMAP2 -ax map-hifi -N 20 -p 0.5 --secondary=yes \
        -t "$THREADS" \
        "$fasta" "${reads[@]}" \
        2> "$ALN_DIR/${key}.aln.log" \
    | $SAMTOOLS sort -@ "$THREADS" -o "$bam"
    $SAMTOOLS index "$bam"
    echo "[align] done: $bam"
}

align_one "$ASM1_FASTA" "${LABEL}.${ASM1_NAME}" &
pid1=$!
align_one "$ASM2_FASTA" "${LABEL}.${ASM2_NAME}" &
pid2=$!
wait $pid1 || { echo "ERROR: alignment failed for ${LABEL}.${ASM1_NAME}"; exit 1; }
wait $pid2 || { echo "ERROR: alignment failed for ${LABEL}.${ASM2_NAME}"; exit 1; }

# ── step 2: bgzip + index each FASTA (skip if bgz + fai exist) ───────────────
bgzip_fasta() {
    local fasta="$1"
    local out="$2"

    if [[ -f "$out" && -f "${out}.fai" ]]; then
        echo "[bgzip] skip $(basename "$fasta"): $out already exists"
        return
    fi

    echo "[bgzip] $(basename "$fasta") → $out"
    $BGZIP -c --threads "$THREADS" "$fasta" > "$out"
    $SAMTOOLS faidx "$out"
    echo "[bgzip] indexed: $out"
}

bgzip_fasta "$ASM1_FASTA" "$BGZ_DIR/${LABEL}.${ASM1_NAME}.fasta.bgz"
bgzip_fasta "$ASM2_FASTA" "$BGZ_DIR/${LABEL}.${ASM2_NAME}.fasta.bgz"

# ── step 3: build SQLite alignment DB (skip assemblies already present) ───────
DB="$DB_DIR/${LABEL}.aln.sqlite"

db_has_asm() {
    local asm="$1"
    [[ -f "$DB" ]] && \
        [[ "$($PYTHON -c "
import sqlite3, sys
try:
    con = sqlite3.connect('$DB')
    n = con.execute('SELECT COUNT(*) FROM aln WHERE assembly=?', ('$asm',)).fetchone()[0]
    print(n)
    con.close()
except Exception:
    print(0)
")" -gt 0 ]]
}

load_asm() {
    local bam="$1"
    local asm="$2"
    if db_has_asm "$asm"; then
        echo "[db] skip ${asm}: already in $DB"
    else
        echo "[db] loading ${asm} from $(basename "$bam")"
        $PYTHON -m qatool.build_aln_db "$bam" "$DB" "$asm"
    fi
}

load_asm "$ALN_DIR/${LABEL}.${ASM1_NAME}.sorted.bam" "$ASM1_NAME"
load_asm "$ALN_DIR/${LABEL}.${ASM2_NAME}.sorted.bam" "$ASM2_NAME"
echo "[db] done: $DB"

echo ""
echo "=== Done: $LABEL ==="
echo "  BAMs:   $ALN_DIR/${LABEL}.{${ASM1_NAME},${ASM2_NAME}}.sorted.bam"
echo "  FASTAs: $BGZ_DIR/${LABEL}.{${ASM1_NAME},${ASM2_NAME}}.fasta.bgz"
echo "  DB:     $DB"
