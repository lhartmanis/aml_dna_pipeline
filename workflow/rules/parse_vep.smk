import pandas as pd


def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


def get_all_samples_from_manifest():
    manifest = pd.read_csv(manifest_path_from_checkpoint(), sep="\t")
    return sorted(manifest["sample_id"].unique().tolist())


rule parse_vep_vcfs:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        vcf_gz=lambda wc: expand(
            "results/vep/{sample}.filtered.vep.vcf.gz",
            sample=get_all_samples_from_manifest()
        ),
        vcf_tbi=lambda wc: expand(
            "results/vep/{sample}.filtered.vep.vcf.gz.tbi",
            sample=get_all_samples_from_manifest()
        ),
        panel_size="results/panel_metrics/panel_size.txt"
    output:
        long="results/analysis/variants_long_annotated.tsv.gz",
        summary="results/analysis/sample_mutect2_summary_annotated.tsv"
    log:
        "logs/snakemake/parse_vep.log"
    shell:
        r"""
        mkdir -p results/analysis logs/snakemake
        python workflow/scripts/11_parse_vep_vcfs.py \
            --vep-dir results/vep \
            --panel-size-file {input.panel_size} \
            --output-long {output.long} \
            --output-summary {output.summary} \
            > {log} 2>&1
        """
