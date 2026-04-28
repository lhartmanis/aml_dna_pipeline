checkpoint build_fastq_manifest:
    output:
        manifest="results/manifest/fastq_manifest.tsv",
        problems="results/manifest/fastq_manifest_problems.tsv",
        summary="results/manifest/fastq_manifest_summary.txt"
    log:
        "logs/build_fastq_manifest.log"
    shell:
        r"""
        mkdir -p results/manifest logs
        python workflow/scripts/01_build_fastq_manifest.py \
            --input-dir {config[input][raw_fastq_dir]} \
            --manifest-out {output.manifest} \
            --problems-out {output.problems} \
            --summary-out {output.summary} \
            > {log} 2>&1
        """
