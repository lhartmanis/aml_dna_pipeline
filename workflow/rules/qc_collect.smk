def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


def get_all_samples_from_checkpoint():
    import pandas as pd
    manifest = pd.read_csv(manifest_path_from_checkpoint(), sep="\t")
    return sorted(manifest["sample_id"].unique().tolist())


rule collect_bam_qc_summary:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        summaries=lambda wc: expand(
            "results/qc/{sample}.qc_summary.json",
            sample=get_all_samples_from_checkpoint()
        )
    output:
        "results/qc/bam_qc_summary.tsv"
    log:
        "logs/snakemake/collect_bam_qc_summary.log"
    shell:
        r"""
        mkdir -p results/qc logs/snakemake
        python workflow/scripts/06_collect_bam_qc_summary.py \
            --qc-dir results/qc \
            --output results/qc/bam_qc_summary.tsv \
            > {log} 2>&1
        """
