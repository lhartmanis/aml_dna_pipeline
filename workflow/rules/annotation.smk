def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


rule annotate_with_vep:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        vcf="results/mutect2_norm/{sample}.filtered.norm.vcf.gz",
        tbi="results/mutect2_norm/{sample}.filtered.norm.vcf.gz.tbi"
    output:
        vcf="results/vep/{sample}.filtered.vep.vcf.gz",
        tbi="results/vep/{sample}.filtered.vep.vcf.gz.tbi"
    log:
        "logs/snakemake/vep_{sample}.log"
    shell:
        r"""
        mkdir -p results/vep logs/vep logs/snakemake
        python workflow/scripts/10_annotate_with_vep.py \
            --input-dir results/mutect2_norm \
            --output-dir results/vep \
            --log-dir logs/vep \
            --sample-id {wildcards.sample} \
            --vep-cache {config[resources][vep_cache]} \
            --vep-fasta {config[resources][vep_fasta]} \
            --assembly GRCh38 \
            --vep-exe {config[tools][vep]} \
            > {log} 2>&1
        """
