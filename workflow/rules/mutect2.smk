def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


rule call_and_filter_mutect2:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        bam="results/bqsr/{sample}.bqsr.bam",
        bai="results/bqsr/{sample}.bqsr.bam.bai"
    output:
        unfiltered_vcf="results/mutect2/{sample}.unfiltered.vcf.gz",
        unfiltered_tbi="results/mutect2/{sample}.unfiltered.vcf.gz.tbi",
        filtered_vcf="results/mutect2/{sample}.filtered.vcf.gz",
        filtered_tbi="results/mutect2/{sample}.filtered.vcf.gz.tbi",
        stats="results/mutect2/{sample}.unfiltered.vcf.gz.stats"
    log:
        "logs/snakemake/mutect2_{sample}.log"
    shell:
        r"""
        mkdir -p results/mutect2 logs/mutect2 logs/snakemake
        python workflow/scripts/08_call_and_filter_mutect2_one_sample.py \
            --bqsr-dir results/bqsr \
            --mutect-dir results/mutect2 \
            --log-dir logs/mutect2 \
            --reference {config[reference][fasta]} \
            --gnomad {config[resources][germline_resource]} \
            --pon {config[resources][pon]} \
            --intervals {config[reference][interval_list]} \
            --sample-id {wildcards.sample} \
            --java-options={config[params][gatk_java_opts]} \
            > {log} 2>&1
        """
