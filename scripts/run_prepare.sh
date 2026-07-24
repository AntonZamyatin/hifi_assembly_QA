#!/bin/bash
set -euo pipefail

# Wrapper around prepare_data.sh — processes multiple species in parallel.
#
# TSV format (tab-separated, no header):
#   label <TAB> reads_dir
#
# Ref assembly: concatenation of the two haplotype FASTAs from examples/genomes/
# LJA assembly: examples/genomes/{label}.lja.fasta

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
EXAMPLES_DIR="$REPO_ROOT/examples"
GENOMES_DIR="$EXAMPLES_DIR/genomes"
ALN_DIR="$EXAMPLES_DIR/alignments"
DB_DIR="$EXAMPLES_DIR/dbs"
BGZ_DIR="$EXAMPLES_DIR/bgz"
TMP_DIR="$EXAMPLES_DIR/tmp"

JOBS=1
THREADS=40
LABELS_FILE=""
TSV_FILE=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --tsv <label_reads.tsv> [--labels <labels.txt>] \\
                         [--jobs <n>] [--threads <n>]

--tsv      Tab-separated file: label<TAB>reads_dir (no header)
--labels   Optional file with one label per line (subset of TSV to process)
           If omitted, all labels in the TSV are processed.
--jobs     Number of labels to process in parallel (default: 1)
--threads  Threads per prepare_data.sh run (default: 40)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tsv)     TSV_FILE="$2";    shift 2 ;;
        --labels)  LABELS_FILE="$2"; shift 2 ;;
        --jobs)    JOBS="$2";        shift 2 ;;
        --threads) THREADS="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

[[ -z "${TSV_FILE:-}" ]] && usage
[[ ! -f "$TSV_FILE" ]] && { echo "ERROR: TSV not found: $TSV_FILE"; exit 1; }

mkdir -p "$ALN_DIR" "$DB_DIR" "$BGZ_DIR" "$TMP_DIR"

# ── load TSV into associative array: label → reads_dir ───────────────────────
declare -A READS_MAP
while IFS=$'\t' read -r label reads_dir; do
    [[ -z "$label" || "$label" == \#* ]] && continue
    READS_MAP["$label"]="$reads_dir"
done < "$TSV_FILE"

# ── determine which labels to process ────────────────────────────────────────
if [[ -n "$LABELS_FILE" ]]; then
    mapfile -t LABELS < "$LABELS_FILE"
else
    LABELS=("${!READS_MAP[@]}")
fi

# ── find haplotype FASTAs for a label (tries common tag conventions) ──────────
find_hap_fastas() {
    local label="$1"
    local hap1="" hap2=""
    for tag1 in hap1 pri mat; do
        local cand="$GENOMES_DIR/${label}.${tag1}.fasta.gz"
        if [[ -f "$cand" ]]; then
            hap1="$cand"
            break
        fi
    done
    for tag2 in hap2 alt pat; do
        local cand="$GENOMES_DIR/${label}.${tag2}.fasta.gz"
        if [[ -f "$cand" ]]; then
            hap2="$cand"
            break
        fi
    done
    echo "$hap1" "$hap2"
}

# ── check if a label is already fully done ───────────────────────────────────
is_done() {
    local label="$1"
    local db="$DB_DIR/${label}.aln.sqlite"
    [[ ! -f "$db" ]] && return 1
    python3 - "$db" <<'PYEOF'
import sqlite3, sys
db = sys.argv[1]
try:
    con = sqlite3.connect(db)
    asms = {r[0] for r in con.execute("SELECT DISTINCT assembly FROM aln")}
    con.close()
    sys.exit(0 if {"ref", "lja"} <= asms else 1)
except Exception:
    sys.exit(1)
PYEOF
}

# ── timestamp helper ──────────────────────────────────────────────────────────
ts() { date '+%H:%M:%S'; }

# ── process one label ─────────────────────────────────────────────────────────
process_label() {
    local label="$1"
    local log="$TMP_DIR/${label}.log"

    echo "[$(ts)] [${label}] starting"

    if is_done "$label"; then
        echo "[$(ts)] [${label}] SKIP — already fully done (DB has ref+lja)"
        return
    fi

    # reads dir
    local reads_dir="${READS_MAP[$label]:-}"
    if [[ -z "$reads_dir" ]]; then
        echo "[$(ts)] [${label}] ERROR: no reads_dir in TSV, skipping" >&2
        return
    fi
    if [[ ! -e "$reads_dir" ]]; then
        echo "[$(ts)] [${label}] ERROR: reads_dir not found: $reads_dir, skipping" >&2
        return
    fi

    # ref haplotype FASTAs
    read -r hap1 hap2 <<< "$(find_hap_fastas "$label")"
    if [[ -z "$hap1" || -z "$hap2" ]]; then
        echo "[$(ts)] [${label}] ERROR: could not find both haplotype FASTAs in $GENOMES_DIR, skipping" >&2
        return
    fi

    # LJA assembly
    local lja="$GENOMES_DIR/${label}.lja.fasta"
    if [[ ! -f "$lja" ]]; then
        echo "[$(ts)] [${label}] ERROR: LJA assembly not found: $lja, skipping" >&2
        return
    fi

    # concatenate haplotypes into diploid ref (skip if already done)
    local ref_fasta="$TMP_DIR/${label}.ref.fasta"
    if [[ ! -f "$ref_fasta" ]]; then
        echo "[$(ts)] [${label}] concat ref: $(basename "$hap1") + $(basename "$hap2")"
        zcat "$hap1" "$hap2" > "$ref_fasta"
        echo "[$(ts)] [${label}] concat ref done"
    else
        echo "[$(ts)] [${label}] SKIP concat ref — already exists"
    fi

    # run prepare_data.sh — tee output to label log AND stdout (prefixed)
    echo "[$(ts)] [${label}] running prepare_data.sh"
    bash "$SCRIPT_DIR/prepare_data.sh" \
        --asm1      "$ref_fasta"  --asm1-name ref \
        --asm2      "$lja"        --asm2-name lja \
        --reads     "$reads_dir" \
        --label     "$label" \
        --aln-dir   "$ALN_DIR" \
        --db-dir    "$DB_DIR" \
        --bgz-dir   "$BGZ_DIR" \
        --threads   "$THREADS" \
        2>&1 | tee "$log" | sed -u "s/^/[${label}] /"

    echo "[$(ts)] [${label}] DONE"
}

# ── run labels with a job pool ────────────────────────────────────────────────
running=0
failed=()

for label in "${LABELS[@]}"; do
    process_label "$label" &
    (( ++running ))
    if (( running >= JOBS )); then
        wait -n || true
        (( --running ))
    fi
done

wait

echo ""
echo "[$(ts)] === All labels processed ==="
