#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 <sample_id_list.txt> [java_opts] [extra_args]"
  echo 'Example: $0 metadata/sample_chunks/chunk_00 -Xmx4g'
  exit 1
fi

SAMPLE_LIST="$1"
JAVA_OPTS="${2:--Xmx4g}"
EXTRA_ARGS="${3:-}"

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
SCRIPT="${PROJECT_DIR}/scripts/mark_duplicates.py"
LOG_DIR="${PROJECT_DIR}/results/logs"

mkdir -p "${LOG_DIR}"

CHUNK_NAME="$(basename "${SAMPLE_LIST}")"
MASTER_LOG="${LOG_DIR}/${CHUNK_NAME}.markdup.master.log"
DONE_TSV="${LOG_DIR}/${CHUNK_NAME}.markdup.done.tsv"
FAIL_TSV="${LOG_DIR}/${CHUNK_NAME}.markdup.failed.tsv"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

echo "=== $(timestamp) START ${CHUNK_NAME} ===" | tee -a "${MASTER_LOG}"

while IFS= read -r sample_id || [[ -n "${sample_id}" ]]; do
  [[ -z "${sample_id}" ]] && continue

  echo "[$(timestamp)] START ${sample_id}" | tee -a "${MASTER_LOG}"

  if python "${SCRIPT}" \
      --sample-id "${sample_id}" \
      --java-options="${JAVA_OPTS}" \
      ${EXTRA_ARGS} >> "${MASTER_LOG}" 2>&1; then
    echo "[$(timestamp)] DONE ${sample_id}" | tee -a "${MASTER_LOG}"
    echo -e "$(timestamp)\t${sample_id}" >> "${DONE_TSV}"
  else
    echo "[$(timestamp)] FAIL ${sample_id}" | tee -a "${MASTER_LOG}"
    echo -e "$(timestamp)\t${sample_id}" >> "${FAIL_TSV}"
  fi

done < "${SAMPLE_LIST}"

echo "=== $(timestamp) END ${CHUNK_NAME} ===" | tee -a "${MASTER_LOG}"
