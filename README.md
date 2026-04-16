# AML DNA panel pipeline

A Snakemake-based workflow for processing targeted DNA sequencing data from AML patient samples, starting from FASTQ files and ending with filtered and annotated somatic variant calls plus sample-level summary tables.

## Overview

This pipeline is designed for targeted AML DNA panel sequencing and currently implements the following major steps:

1. Build a FASTQ manifest from raw sequencing files
2. Align reads per sequencing unit with bwa-mem2
3. Merge unit-level BAMs per sample
4. Mark duplicates
5. Compute BAM QC and coverage summaries
6. Perform base quality score recalibration (BQSR)
7. Call somatic variants with GATK Mutect2
8. Filter Mutect2 calls
9. Normalize filtered VCFs
10. Annotate variants with Ensembl VEP
11. Parse PASS variants into tabular sample-level summaries
12. Parse VEP-annotated variants into analysis-ready tables

## Workflow structure

- `workflow/Snakefile` – main Snakemake entry point
- `workflow/rules/` – rule files split by workflow stage
- `workflow/scripts/` – Python and shell scripts used by rules
- `config/` – configuration files, sample sheet template, panel files
- `resources/` – notes on references and downloads
- `profiles/cluster/` – cluster execution profile
- `results/` – output directory
- `logs/` – log directory
- `benchmark/` – benchmark/runtime files
- `dev/` – helper scripts and legacy launcher scripts retained for development/debugging

## Inputs

The workflow expects:

- paired-end FASTQ files
- a sample sheet describing samples and FASTQ locations
- panel target intervals in BED format
- panel gene list
- reference genome and associated indices/resources
- external resources for BQSR and somatic calling

## Main outputs

Typical outputs include:

- sample FASTQ manifest
- aligned, merged, duplicate-marked, and recalibrated BAM files
- BAM QC summary tables
- filtered and normalized VCF files
- VEP-annotated VCF files
- parsed variant summary tables for downstream analysis

## Software

Main tools used in this workflow include:

- Snakemake
- bwa-mem2
- samtools
- Picard
- GATK
- bcftools
- Ensembl VEP
- Python

Tool versions should be recorded in the manuscript and/or release notes for each analysis run.

## Configuration

Main runtime settings are stored in:

- `config/config.yaml`

Panel resources are stored in:

- `config/panel/panel_targets.bed`
- `config/panel/panel_genes.txt`

A template sample sheet is provided in:

- `config/samples.example.tsv`

## Running the workflow

### Dry run

```snakemake -s workflow/Snakefile -n```

### Local execution

```snakemake -s workflow/Snakefile --cores 8```

### Cluster/profile execution

```snakemake -s workflow/Snakefile --profile profiles/cluster```

### Notes

- Large genomics files and generated outputs are excluded from version control via .gitignore.
- Empty output directories are retained using .gitkeep.
- Development/helper scripts from the original project are stored in dev/ and are not required for normal workflow execution.

### Status

This repository is under active development as the original project scripts are being refactored into a reproducible Snakemake workflow.
