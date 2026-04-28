#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path


DEFAULT_SAMPLE_BAM_DIR = Path("results/markdup")
DEFAULT_BQSR_DIR = Path("results/bqsr")
DEFAULT_LOG_DIR = Path("logs/bqsr")
DEFAULT_REF_FA = Path("resources/reference.fa")
DEFAULT_DBSNP = Path("resources/dbsnp.vcf.gz")
DEFAULT_MILLS = Path("resources/mills.vcf.gz")


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfh) -> None:
    line = f"[{timestamp()}] {msg}"
    print(line)
    logfh.write(line + "\n")
    logfh.flush()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run GATK BQSR on one sample.")
    parser.add_argument("--sample-bam-dir", default=str(DEFAULT_SAMPLE_BAM_DIR))
    parser.add_argument("--bqsr-dir", default=str(DEFAULT_BQSR_DIR))
    parser.add_argument("--log-dir", default=str(DEFAULT_LOG_DIR))
    parser.add_argument("--reference", default=str(DEFAULT_REF_FA))
    parser.add_argument("--dbsnp", default=str(DEFAULT_DBSNP))
    parser.add_argument("--mills", default=str(DEFAULT_MILLS))
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--plot", action="store_true", help="Run AnalyzeCovariates and make PDF plots")
    parser.add_argument("--java-options", default="-Xmx8g")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    sample_bam_dir = Path(args.sample_bam_dir)
    bqsr_dir = Path(args.bqsr_dir)
    log_dir = Path(args.log_dir)
    reference = Path(args.reference)
    dbsnp = Path(args.dbsnp)
    mills = Path(args.mills)

    for p in [sample_bam_dir, reference, dbsnp, mills]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required file/path: {p}")

    bqsr_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    sample_id = args.sample_id

    input_bam = sample_bam_dir / f"{sample_id}.markdup.bam"
    input_bai = sample_bam_dir / f"{sample_id}.markdup.bam.bai"

    # tolerate either Picard- or samtools-style index naming
    if not input_bai.exists():
        alt_bai = sample_bam_dir / f"{sample_id}.markdup.bai"
        if alt_bai.exists():
            input_bai = alt_bai

    pre_table = bqsr_dir / f"{sample_id}.pre_recal.table"
    output_bam = bqsr_dir / f"{sample_id}.bqsr.bam"
    output_bai = bqsr_dir / f"{sample_id}.bqsr.bam.bai"
    post_table = bqsr_dir / f"{sample_id}.post_recal.table"
    plots_pdf = bqsr_dir / f"{sample_id}.bqsr_plots.pdf"

    sample_log = log_dir / f"{sample_id}.bqsr.log"

    if not input_bam.exists():
        raise FileNotFoundError(f"Missing input BAM: {input_bam}")
    if not input_bai.exists():
        raise FileNotFoundError(f"Missing input BAI: {input_bai}")

    expected_outputs = [pre_table, output_bam, post_table]
    if args.plot:
        expected_outputs.append(plots_pdf)

    if all(p.exists() for p in expected_outputs) and not args.force:
        with open(sample_log, "a") as logfh:
            log(f"SKIP {sample_id}: outputs already exist", logfh)
        return

    base_recal_cmd = [
        "gatk", "--java-options", args.java_options,
        "BaseRecalibrator",
        "-R", str(reference),
        "-I", str(input_bam),
        "--known-sites", str(dbsnp),
        "--known-sites", str(mills),
        "-O", str(pre_table),
    ]

    apply_bqsr_cmd = [
        "gatk", "--java-options", args.java_options,
        "ApplyBQSR",
        "-R", str(reference),
        "-I", str(input_bam),
        "--bqsr-recal-file", str(pre_table),
        "-O", str(output_bam),
    ]

    index_cmd = ["samtools", "index", str(output_bam)]

    post_recal_cmd = [
        "gatk", "--java-options", args.java_options,
        "BaseRecalibrator",
        "-R", str(reference),
        "-I", str(output_bam),
        "--known-sites", str(dbsnp),
        "--known-sites", str(mills),
        "-O", str(post_table),
    ]

    analyze_cmd = [
        "gatk", "--java-options", args.java_options,
        "AnalyzeCovariates",
        "-before", str(pre_table),
        "-after", str(post_table),
        "-plots", str(plots_pdf),
    ]

    cmds = [
        ("BaseRecalibrator(pre)", base_recal_cmd),
        ("ApplyBQSR", apply_bqsr_cmd),
        ("samtools index", index_cmd),
        ("BaseRecalibrator(post)", post_recal_cmd),
    ]
    if args.plot:
        cmds.append(("AnalyzeCovariates", analyze_cmd))

    if args.dry_run:
        for label, cmd in cmds:
            print(f"\n# {label}")
            print(" ".join(shlex.quote(x) for x in cmd))
        return

    start = time.time()

    with open(sample_log, "a") as logfh:
        logfh.write(f"Start time: {datetime.now()}\n")
        logfh.write(f"Sample: {sample_id}\n")
        logfh.write(f"Input BAM: {input_bam}\n")
        logfh.write(f"Reference: {reference}\n")
        logfh.write(f"dbSNP: {dbsnp}\n")
        logfh.write(f"Mills: {mills}\n")
        logfh.write(f"Plot mode: {args.plot}\n\n")
        logfh.flush()

        try:
            for label, cmd in cmds:
                log(f"RUN {label}", logfh)
                logfh.write("Command:\n")
                logfh.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
                logfh.flush()
                subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh)

            # final sanity checks
            for p in [pre_table, output_bam, post_table]:
                if not p.exists():
                    raise RuntimeError(f"Expected output missing: {p}")

            if not output_bai.exists():
                raise RuntimeError(f"Expected BAM index missing: {output_bai}")

            if args.plot and not plots_pdf.exists():
                raise RuntimeError(f"Expected plot PDF missing: {plots_pdf}")

            elapsed = time.time() - start
            logfh.write(f"\nEnd time: {datetime.now()}\n")
            logfh.write(f"Elapsed seconds: {elapsed:.2f}\n")
            logfh.flush()
            log(f"DONE {sample_id}: elapsed={elapsed:.2f}s", logfh)

        except Exception as e:
            log(f"FAIL {sample_id}: {e}", logfh)
            raise


if __name__ == "__main__":
    main()
