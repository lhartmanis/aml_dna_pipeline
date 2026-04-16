#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/leonard.hartmanis/proj/DNA_seq"
INPUT_DIR="${PROJECT_DIR}/results/mutect2_norm"
OUT_DIR="${PROJECT_DIR}/results/vep"
META_DIR="${PROJECT_DIR}/metadata"

SAMPLE_LIST="${META_DIR}/all_sample_ids.txt"

VEP_DIR="/home/leonard.hartmanis/packages/ensembl-vep"
VEP_CACHE_DIR="/home/leonard.hartmanis/.vep"
VEP_FASTA="${VEP_CACHE_DIR}/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

mkdir -p "${OUT_DIR}"

# Needed for Bio::DB::HTS to find libhts.so.1 in your conda env
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list not found: ${SAMPLE_LIST}" >&2
    exit 1
fi

if [[ ! -x "${VEP_DIR}/vep" ]]; then
    echo "ERROR: VEP executable not found: ${VEP_DIR}/vep" >&2
    exit 1
fi

if [[ ! -f "${VEP_FASTA}" ]]; then
    echo "ERROR: VEP FASTA not found: ${VEP_FASTA}" >&2
    exit 1
fi

while IFS= read -r sample_id; do
    [[ -z "${sample_id}" ]] && continue

    in_vcf="${INPUT_DIR}/${sample_id}.filtered.norm.vcf.gz"
    out_vcf="${OUT_DIR}/${sample_id}.filtered.vep.vcf.gz"

    if [[ ! -f "${in_vcf}" ]]; then
        echo "WARNING: missing input VCF for ${sample_id}: ${in_vcf}" >&2
        continue
    fi

    if [[ -f "${out_vcf}" && -f "${out_vcf}.tbi" ]]; then
        echo "Skipping ${sample_id} (already annotated)"
        continue
    fi

    echo "Annotating ${sample_id}..."

    "${VEP_DIR}/vep" \
        --cache \
        --offline \
        --assembly GRCh38 \
        --dir_cache "${VEP_CACHE_DIR}" \
        --fasta "${VEP_FASTA}" \
        --input_file "${in_vcf}" \
        --output_file "${out_vcf}" \
        --vcf \
        --compress_output bgzip \
        --force_overwrite \
        --symbol \
        --hgvs \
        --terms SO \
        --pick \
        --pick_order canonical,appris,biotype,rank,length \
        --everything

    tabix -f -p vcf "${out_vcf}"
done < "${SAMPLE_LIST}"

echo "All done."
