#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <rgid_list.txt> [bwa_threads] [sort_threads]"
    exit 1
fi

LIST="$1"
BWA_THREADS="${2:-6}"
SORT_THREADS="${3:-3}"

PROJECT="/home/leonard.hartmanis/proj/DNA_seq"
SCRIPT="${PROJECT}/scripts/align_one_unit.py"
UNIT_BAM_DIR="${PROJECT}/results/unit_bam"
LOG_DIR="${PROJECT}/results/logs"
RUN_NAME="$(basename "$LIST")"

MASTER_LOG="${LOG_DIR}/${RUN_NAME}.master.log"
FAIL_LOG="${LOG_DIR}/${RUN_NAME}.failures.tsv"
DONE_LOG="${LOG_DIR}/${RUN_NAME}.done.tsv"

mkdir -p "$LOG_DIR"

echo "=== $(date) START ${RUN_NAME} ===" | tee -a "$MASTER_LOG"
echo -e "timestamp\trgid\treason" >> "$FAIL_LOG"
echo -e "timestamp\trgid\tbam" >> "$DONE_LOG"

while IFS= read -r RGID || [ -n "$RGID" ]; do
    [ -z "$RGID" ] && continue

    BAM="${UNIT_BAM_DIR}/${RGID}.sorted.bam"
    BAI="${UNIT_BAM_DIR}/${RGID}.sorted.bam.bai"

    if [ -s "$BAM" ] && [ -s "$BAI" ]; then
        echo "[$(date)] SKIP ${RGID} (already done)" | tee -a "$MASTER_LOG"
        continue
    fi

    echo "[$(date)] START ${RGID}" | tee -a "$MASTER_LOG"

    if python "$SCRIPT" \
        --rgid "$RGID" \
        --bwa-threads "$BWA_THREADS" \
        --sort-threads "$SORT_THREADS"
    then
        echo "[$(date)] DONE ${RGID}" | tee -a "$MASTER_LOG"
        echo -e "$(date +'%F %T')\t${RGID}\t${BAM}" >> "$DONE_LOG"
    else
        echo "[$(date)] FAIL ${RGID}" | tee -a "$MASTER_LOG"
        echo -e "$(date +'%F %T')\t${RGID}\talign_one_unit_failed" >> "$FAIL_LOG"
    fi
done < "$LIST"

echo "=== $(date) END ${RUN_NAME} ===" | tee -a "$MASTER_LOG"
