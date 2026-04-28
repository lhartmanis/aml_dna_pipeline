def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


rule bam_qc_summary:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        bam="results/markdup/{sample}.markdup.bam",
        bai="results/markdup/{sample}.markdup.bam.bai"
    output:
        flagstat="results/qc/{sample}.flagstat.txt",
        idxstats="results/qc/{sample}.idxstats.txt",
        depth="results/qc/{sample}.target_depth.tsv",
        hsmetrics="results/qc/{sample}.hs_metrics.txt",
        summary="results/qc/{sample}.qc_summary.json"
    log:
        "logs/snakemake/qc_{sample}.log"
    shell:
        r"""
        mkdir -p results/qc logs/qc logs/snakemake
        python workflow/scripts/05_bam_qc_summary.py \
            --sample-bam-dir results/markdup \
            --qc-dir results/qc \
            --log-dir logs/qc \
            --interval-bed {config[reference][interval_bed]} \
            --interval-list {config[reference][interval_list]} \
            --reference {config[reference][fasta]} \
            --sample-id {wildcards.sample} \
            > {log} 2>&1
        """
