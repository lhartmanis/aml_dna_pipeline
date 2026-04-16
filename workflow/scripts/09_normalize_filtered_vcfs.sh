#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
RESULTS_DIR="${PROJECT_DIR}/results"
REF_DIR="${PROJECT_DIR}/ref"
META_DIR="${PROJECT_DIR}/metadata"

# Adjust this if your fasta has a different name
REF_FASTA="${REF_DIR}/GRCh38.p14.genome.fa"

INPUT_DIR="${RESULTS_DIR}/mutect2"
OUT_DIR="${RESULTS_DIR}/mutect2_norm"
SAMPLE_LIST="${META_DIR}/all_sample_ids.txt"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${REF_FASTA}" ]]; then
    echo "ERROR: Reference fasta not found: ${REF_FASTA}" >&2
    exit 1
fi

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}" >&2
    exit 1
fi

while IFS= read -r sample_id; do
    [[ -z "${sample_id}" ]] && continue

    in_vcf="${INPUT_DIR}/${sample_id}.filtered.vcf.gz"
    out_vcf="${OUT_DIR}/${sample_id}.filtered.norm.vcf.gz"

    if [[ ! -f "${in_vcf}" ]]; then
        echo "WARNING: Missing input VCF for ${sample_id}: ${in_vcf}" >&2
        continue
    fi

    echo "Normalizing ${sample_id}..."

    bcftools norm \
        -m -any \
        -f "${REF_FASTA}" \
        -Oz \
        -o "${out_vcf}" \
        "${in_vcf}"

    tabix -f -p vcf "${out_vcf}"
done < "${SAMPLE_LIST}"

echo "Done."
