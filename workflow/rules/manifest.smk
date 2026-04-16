from pathlib import Path

REPO_ROOT = Path(workflow.basedir).parent

rule build_fastq_manifest:
    output:
        manifest = str(REPO_ROOT / "results" / "manifest" / "fastq_manifest.tsv"),
        problems = str(REPO_ROOT / "results" / "manifest" / "fastq_manifest_problems.tsv"),
        summary = str(REPO_ROOT / "results" / "manifest" / "fastq_manifest_summary.txt"),
    params:
        script = str(REPO_ROOT / "workflow" / "scripts" / "01_build_fastq_manifest.py"),
        input_dir = config["input"]["raw_fastq_dir"],
    log:
        str(REPO_ROOT / "logs" / "build_fastq_manifest.log")
    conda:
        str(REPO_ROOT / "workflow" / "envs" / "dna_panel.yaml")
    shell:
        r"""
        mkdir -p {REPO_ROOT}/results/manifest {REPO_ROOT}/logs
        python {params.script} \
          --input-dir {params.input_dir} \
          --manifest-out {output.manifest} \
          --problems-out {output.problems} \
          --summary-out {output.summary} \
          > {log} 2>&1
        """
