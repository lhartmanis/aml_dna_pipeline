#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
INTERVALS_DIR="${PROJECT_DIR}/intervals"
OUT_DIR="${PROJECT_DIR}/results/panel_metrics"
mkdir -p "${OUT_DIR}"

BED_FILE="${INTERVALS_DIR}/panel_targets_all_hg38.merged.bed"
OUT_TXT="${OUT_DIR}/panel_size.txt"

if [[ ! -f "${BED_FILE}" ]]; then
    echo "ERROR: BED file not found: ${BED_FILE}" >&2
    exit 1
fi

awk 'BEGIN{sum=0} {sum += ($3-$2)} END {printf "panel_bp\t%d\npanel_mb\t%.6f\n", sum, sum/1e6}' \
    "${BED_FILE}" > "${OUT_TXT}"

cat "${OUT_TXT}"
