# 🧬 AML DNA Panel Pipeline 

A modular Snakemake workflow for processing targeted DNA sequencing data from acute myeloid leukemia (AML) samples, from raw FASTQ files to filtered, annotated somatic variant calls and downstream analysis-ready summary tables.

## Quick start

```bash
git clone https://github.com/lhartmanis/aml_dna_panel_pipeline.git
cd aml_dna_panel_pipeline

conda env create -f workflow/envs/dna_panel.yaml
conda activate dna_panel

# preview the full workflow without running jobs
snakemake -s workflow/Snakefile -np --cores 1
```

## What the pipeline does

The workflow implements the following stages:

1. Build a FASTQ manifest from sample directories
2. Align reads per sequencing unit with `bwa-mem2`
3. Merge lane/unit BAMs per sample
4. Mark duplicates with Picard
5. Generate BAM QC and coverage summaries
6. Perform base quality score recalibration (BQSR) with GATK
7. Call somatic variants with GATK Mutect2 and filter calls
8. Normalize filtered VCFs with `bcftools`
9. Annotate variants with Ensembl VEP
10. Parse annotated variants into analysis-ready long and summary tables

## Repository layout

```text
aml_dna_panel_pipeline/
├── config/                 # user-editable configuration files
├── dev/                    # helper and development scripts
├── profiles/               # optional execution profiles
├── resources/              # notes on references and resource setup
├── workflow/               # Snakefile, rule files, scripts, environments
├── LICENSE
└── README.md
```

## Requirements

- Conda or Mamba
- Linux or macOS shell environment
- Local copies of the reference genome and somatic-calling resources
- Sufficient disk space for BAM/VCF generation

## Installation

Create the main environment:

```bash
conda env create -f workflow/envs/dna_panel.yaml
conda activate dna_panel
```

If you use a separate VEP environment:

```bash
conda env create -f workflow/envs/vep.yaml
```

Check that core tools are available:

```bash
which bwa-mem2
which samtools
which gatk
which picard
which bcftools
snakemake --version
python --version
```

## Reference and resource setup

Before running the workflow, prepare the reference FASTA and required resources locally.

Expected resources include:

- GRCh38 reference FASTA
- FASTA index (`.fai`)
- sequence dictionary (`.dict`)
- dbSNP resource for BQSR
- Mills and 1000G gold-standard indels
- GATK panel of normals
- gnomAD germline resource
- VEP cache
- VEP FASTA for offline annotation

Reference download notes are provided in:

```text
resources/reference_download_notes.md
```

### Reference indexing

Before alignment, the reference FASTA must be indexed for `bwa-mem2`:

```bash
bwa-mem2 index /path/to/GRCh38.p14.genome.fa
```

### VEP cache installation

For offline VEP annotation, install the GRCh38 cache into your VEP cache directory:

```bash
vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to_vep_cache/.vep
```

If your configuration points to a custom VEP cache location, make sure `config/config.yaml` matches that path.

## Configuration

Main runtime settings are stored in:

```text
config/config.yaml
```

Important fields include:

- `input.raw_fastq_dir`: root directory containing per-sample FASTQ folders
- `reference.*`: reference genome and interval resources
- `resources.*`: dbSNP, Mills, PoN, gnomAD, and VEP resources
- `tools.*`: executable names or paths
- `params.*`: thread counts and Java memory settings

Existing `config/config.yaml` paths are placeholders and should be updated to match your local installation.

Panel-specific files are stored in:

```text
config/panel/panel_targets.bed
config/panel/panel_genes.txt
```

A template sample sheet is provided in:

```text
config/samples.example.tsv
```

## Running the workflow

### Dry run the full pipeline

```bash
snakemake -s workflow/Snakefile -np --cores 1
```

### Run the full pipeline locally

```bash
snakemake -s workflow/Snakefile --cores 8
```

### Run with a profile

```bash
snakemake -s workflow/Snakefile --profile profiles/cluster
```

### Run a specific target

For example, build one aligned BAM:

```bash
snakemake -s workflow/Snakefile --cores 1 results/bam/target_sample.sorted.bam
```

Or build the final annotated summary table:

```bash
snakemake -s workflow/Snakefile --cores 8 results/analysis/sample_mutect2_summary_annotated.tsv
```

## Main outputs

The workflow produces the following output categories:

- FASTQ manifest and parsing summaries
- aligned, merged, duplicate-marked, and BQSR-corrected BAMs
- per-sample BAM QC metrics and aggregated QC summary tables
- filtered and normalized Mutect2 VCFs
- VEP-annotated VCFs
- final analysis-ready tables:
  - `results/analysis/variants_long_annotated.tsv.gz`
  - `results/analysis/sample_mutect2_summary_annotated.tsv`

## Workflow notes

- The workflow uses a Snakemake checkpoint to build the FASTQ manifest before sample- and read-group–specific expansion.
- Large intermediate BAM files are marked as temporary where appropriate, so completed downstream stages can trigger cleanup automatically.
- Final downstream variant parsing is driven from VEP-annotated VCFs and a computed panel-size file.
- Large runtime outputs are excluded from version control via `.gitignore`.

## Development status

This repository contains a modular, reproducible implementation of the main AML targeted-panel DNA processing workflow. Helper scripts under `dev/` are retained for development, testing, and legacy debugging, but are not required for standard execution.

## License

This project is released under the MIT License. See `LICENSE`.
