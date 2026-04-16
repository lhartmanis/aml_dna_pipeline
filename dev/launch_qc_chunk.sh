#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
cd "${PROJECT_DIR}"

mkdir -p results/logs   # <-- FIX

for chunk in metadata/sample_chunks_25/chunk_*; do
  [[ -e "$chunk" ]] || { echo "No chunk files found"; exit 1; }
  chunk_name=$(basename "$chunk")
  nohup bash scripts/run_qc_chunk.sh "$chunk" \
    > "results/logs/${chunk_name}.qc.nohup.out" 2>&1 &
done
