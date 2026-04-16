#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
cd "${PROJECT_DIR}"

for chunk in metadata/sample_chunks/chunk_*; do
  [[ -e "$chunk" ]] || { echo "No chunk files found"; exit 1; }
  chunk_name=$(basename "$chunk")
  nohup bash scripts/run_markdup_chunk.sh "$chunk" -Xmx4g \
    > "results/logs/${chunk_name}.markdup.nohup.out" 2>&1 &
done
