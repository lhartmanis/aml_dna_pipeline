#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
REFERENCE=""
SAMPLE_ID=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --sample-id)
            SAMPLE_ID="$2"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "${INPUT_DIR}" || -z "${OUTPUT_DIR}" || -z "${REFERENCE}" || -z "${SAMPLE_ID}" ]]; then
    echo "ERROR: Required arguments are missing." >&2
    echo "Usage: $0 --input-dir DIR --output-dir DIR --reference FASTA --sample-id ID" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

if [[ ! -f "${REFERENCE}" ]]; then
    echo "ERROR: Reference fasta not found: ${REFERENCE}" >&2
    exit 1
fi

in_vcf="${INPUT_DIR}/${SAMPLE_ID}.filtered.vcf.gz"
out_vcf="${OUTPUT_DIR}/${SAMPLE_ID}.filtered.norm.vcf.gz"

if [[ ! -f "${in_vcf}" ]]; then
    echo "ERROR: Missing input VCF for ${SAMPLE_ID}: ${in_vcf}" >&2
    exit 1
fi

echo "Normalizing ${SAMPLE_ID}..."

bcftools norm \
    -m -any \
    -f "${REFERENCE}" \
    -Oz \
    -o "${out_vcf}" \
    "${in_vcf}"

tabix -f -p vcf "${out_vcf}"

echo "Done: ${SAMPLE_ID}"
