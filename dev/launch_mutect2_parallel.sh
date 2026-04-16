#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
cd "${PROJECT_DIR}"

mkdir -p results/logs

for chunk in metadata/sample_chunks_15/chunk_*; do
  [[ -e "$chunk" ]] || { echo "No chunk files found"; exit 1; }
  chunk_name=$(basename "$chunk")
  nohup bash scripts/run_mutect2_chunk.sh "$chunk" -Xmx4g \
    > "results/logs/${chunk_name}.mutect2.nohup.out" 2>&1 &
done
