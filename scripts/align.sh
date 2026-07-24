#!/bin/bash
#set -euo pipefail

REF=$1
READS_INPUT=$2
LABEL=$3

READS=()
if [[ -f "$READS_INPUT" ]]; then
    READS+=("$READS_INPUT")
else
    while IFS= read -r -d '' f; do
        READS+=("$f")
    done < <(find "$READS_INPUT" -maxdepth 1 -type f -name "*.fasta" -print0 | sort -z)
fi

#source /scratch/azz5280/mambaforge/etc/profile.d/conda.sh
#conda activate base
MINIMAP2=/scratch/azz5280/mambaforge/bin/minimap2

#echo $(which minimap2)
echo $(${MINIMAP2} --version)


echo "${READS[@]}"
if [[ ${#READS[@]} -eq 0 ]]; then
    echo "ERROR: No .fasta reads found in $READS_INPUT"
    exit 1
fi

echo "align ${READS[*]} to reference $REF"

echo "=====================\
    minimap2 -ax map-hifi -N 20 -p 0.5 --secondary=yes \
    -t 100 \
    --split-prefix alignments/${LABEL}_index.tmp \
    ${REF} ${READS[*]} \
    2> alignments/${LABEL}.aln.log \
| samtools sort -@ 40 -o alignments/${LABEL}.sorted.bam

samtools index alignments/${LABEL}.sorted.bam\
     ======================"

$MINIMAP2 -ax map-hifi -N 20 -p 0.5 --secondary=yes \
    -t 100 \
    --split-prefix alignments/${LABEL}_index.tmp \
    "${REF}" ${READS[@]} \
    2> alignments/${LABEL}.aln.log \
| samtools sort -@ 40 -o alignments/${LABEL}.sorted.bam

samtools index alignments/${LABEL}.sorted.bam
