def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


rule mark_duplicates:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        bam="results/sample_bam/{sample}.merged.bam",
        bai="results/sample_bam/{sample}.merged.bam.bai"
    output:
        bam=temp("results/markdup/{sample}.markdup.bam"),
        bai=temp("results/markdup/{sample}.markdup.bam.bai"),
        metrics="results/qc/{sample}.markdup.metrics.txt"
    log:
        "logs/snakemake/markdup_{sample}.log"
    shell:
        r"""
        mkdir -p results/markdup results/qc logs/markdup logs/snakemake
        python workflow/scripts/04_mark_duplicates.py \
            --sample-bam-dir results/sample_bam \
            --output-dir results/markdup \
            --qc-dir results/qc \
            --log-dir logs/markdup \
            --sample-id {wildcards.sample} \
            --java-options={config[params][gatk_java_opts]} \
            > {log} 2>&1
        """
