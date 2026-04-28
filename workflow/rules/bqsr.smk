def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


rule bqsr_one_sample:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        bam="results/markdup/{sample}.markdup.bam",
        bai="results/markdup/{sample}.markdup.bam.bai"
    output:
        table_pre="results/bqsr/{sample}.pre_recal.table",
        bam="results/bqsr/{sample}.bqsr.bam",
        bai="results/bqsr/{sample}.bqsr.bam.bai",
        table_post="results/bqsr/{sample}.post_recal.table"
    log:
        "logs/snakemake/bqsr_{sample}.log"
    params:
        ref=config["reference"]["fasta"],
        dbsnp=config["resources"]["dbsnp"],
        mills=config["resources"]["mills_indels"],
        java_opts=config["params"]["gatk_java_opts"]
    shell:
        r"""
        mkdir -p results/bqsr logs/bqsr logs/snakemake
        python workflow/scripts/07_bqsr_one_sample.py \
            --sample-bam-dir results/markdup \
            --bqsr-dir results/bqsr \
            --log-dir logs/bqsr \
            --reference {params.ref} \
            --dbsnp {params.dbsnp} \
            --mills {params.mills} \
            --sample-id {wildcards.sample} \
            --java-options={params.java_opts} \
            > {log} 2>&1
        """
