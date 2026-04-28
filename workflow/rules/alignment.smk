import pandas as pd


def manifest_path_from_checkpoint():
    ckpt = checkpoints.build_fastq_manifest.get()
    return ckpt.output.manifest


def get_row_for_rgid(rgid):
    manifest = pd.read_csv(manifest_path_from_checkpoint(), sep="\t")
    row = manifest.loc[manifest["rgid"] == rgid]
    if row.empty:
        raise ValueError(f"No manifest row found for rgid={rgid}")
    return row.iloc[0]


rule align_one_unit:
    input:
        manifest=lambda wc: manifest_path_from_checkpoint()
    output:
        bam=temp("results/bam/{rgid}.sorted.bam"),
        bai=temp("results/bam/{rgid}.sorted.bam.bai")
    log:
        "logs/snakemake/align_{rgid}.log"
    shell:
        r"""
        mkdir -p results/bam logs/alignment logs/snakemake
        python workflow/scripts/02_align_one_unit.py \
          --rgid {wildcards.rgid} \
          --manifest {input.manifest} \
          --reference {config[reference][fasta]} \
          --outdir results/bam \
          --logdir logs/alignment \
          --bwa-threads {config[params][bwa_threads]} \
          --sort-threads {config[params][samtools_threads]} \
          --sort-memory 1G \
          > {log} 2>&1
        """

