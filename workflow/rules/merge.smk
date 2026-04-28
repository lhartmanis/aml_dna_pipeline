import pandas as pd


def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


def get_rgids_for_sample(sample):
    manifest = pd.read_csv(manifest_path_from_checkpoint(), sep="\t")
    rgids = manifest.loc[manifest["sample_id"] == sample, "rgid"].tolist()
    if not rgids:
        raise ValueError(f"No rgids found for sample={sample}")
    return rgids


rule merge_sample_bams:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint(),
        unit_bams=lambda wc: expand(
            "results/bam/{rgid}.sorted.bam",
            rgid=get_rgids_for_sample(wc.sample)
        ),
        unit_bais=lambda wc: expand(
            "results/bam/{rgid}.sorted.bam.bai",
            rgid=get_rgids_for_sample(wc.sample)
        )
    output:
        bam=temp("results/sample_bam/{sample}.merged.bam"),
        bai=temp("results/sample_bam/{sample}.merged.bam.bai")
    log:
        "logs/snakemake/merge_{sample}.log"
    threads: 4
    shell:
        r"""
        mkdir -p results/sample_bam logs/merge logs/snakemake
        python workflow/scripts/03_merge_sample_bams.py \
            --manifest {input.manifest} \
            --unit-bam-dir results/bam \
            --sample-bam-dir results/sample_bam \
            --log-dir logs/merge \
            --sample-id {wildcards.sample} \
            --samtools-threads {threads} \
            > {log} 2>&1
        """
