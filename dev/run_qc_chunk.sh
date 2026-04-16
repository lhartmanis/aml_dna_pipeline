#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <sample_id_list.txt> [extra_args]"
  echo "Example: $0 /home/leonard.hartmanis/proj/DNA_seq/metadata/sample_chunks/chunk_00"
  echo "Example: $0 /home/leonard.hartmanis/proj/DNA_seq/metadata/sample_chunks/chunk_00 --force"
  exit 1
fi

SAMPLE_LIST="$1"
EXTRA_ARGS="${2:-}"

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
SCRIPT="${PROJECT_DIR}/scripts/bam_qc_summary.py"
LOG_DIR="${PROJECT_DIR}/results/logs"

mkdir -p "${LOG_DIR}"

CHUNK_NAME="$(basename "${SAMPLE_LIST}")"
MASTER_LOG="${LOG_DIR}/${CHUNK_NAME}.qc.master.log"
DONE_TSV="${LOG_DIR}/${CHUNK_NAME}.qc.done.tsv"
FAIL_TSV="${LOG_DIR}/${CHUNK_NAME}.qc.failed.tsv"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

echo "=== $(timestamp) START ${CHUNK_NAME} ===" | tee -a "${MASTER_LOG}"

if [[ ! -f "${SAMPLE_LIST}" ]]; then
  echo "[$(timestamp)] FAIL missing sample list: ${SAMPLE_LIST}" | tee -a "${MASTER_LOG}"
  exit 1
fi

if [[ ! -s "${DONE_TSV}" ]]; then
  echo -e "timestamp\tsample_id" > "${DONE_TSV}"
fi

if [[ ! -s "${FAIL_TSV}" ]]; then
  echo -e "timestamp\tsample_id" > "${FAIL_TSV}"
fi

while IFS= read -r sample_id || [[ -n "${sample_id}" ]]; do
  [[ -z "${sample_id}" ]] && continue

  echo "[$(timestamp)] START ${sample_id}" | tee -a "${MASTER_LOG}"

  if python "${SCRIPT}" \
      --sample-id "${sample_id}" \
      ${EXTRA_ARGS} >> "${MASTER_LOG}" 2>&1; then
    echo "[$(timestamp)] DONE ${sample_id}" | tee -a "${MASTER_LOG}"
    echo -e "$(timestamp)\t${sample_id}" >> "${DONE_TSV}"
  else
    echo "[$(timestamp)] FAIL ${sample_id}" | tee -a "${MASTER_LOG}"
    echo -e "$(timestamp)\t${sample_id}" >> "${FAIL_TSV}"
  fi

done < "${SAMPLE_LIST}"

echo "=== $(timestamp) END ${CHUNK_NAME} ===" | tee -a "${MASTER_LOG}"
