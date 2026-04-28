def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


rule normalize_filtered_vcf:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        vcf="results/mutect2/{sample}.filtered.vcf.gz",
        tbi="results/mutect2/{sample}.filtered.vcf.gz.tbi"
    output:
        vcf="results/mutect2_norm/{sample}.filtered.norm.vcf.gz",
        tbi="results/mutect2_norm/{sample}.filtered.norm.vcf.gz.tbi"
    log:
        "logs/snakemake/norm_{sample}.log"
    shell:
        r"""
        mkdir -p results/mutect2_norm logs/snakemake
        bash workflow/scripts/09_normalize_filtered_vcfs.sh \
            --input-dir results/mutect2 \
            --output-dir results/mutect2_norm \
            --reference {config[reference][fasta]} \
            --sample-id {wildcards.sample} \
            > {log} 2>&1
        """
