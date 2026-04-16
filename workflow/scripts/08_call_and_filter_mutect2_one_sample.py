#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path


DEFAULT_BQSR_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/bqsr")
DEFAULT_MUTECT_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/mutect2")
DEFAULT_LOG_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/logs")

DEFAULT_REF_FA = Path("/home/leonard.hartmanis/proj/DNA_seq/ref/GRCh38.p14.genome.fa")
DEFAULT_GNOMAD = Path("/home/leonard.hartmanis/proj/DNA_seq/ref/af-only-gnomad.hg38.vcf.gz")
DEFAULT_PON = Path("/home/leonard.hartmanis/proj/DNA_seq/ref/1000g_pon.hg38.vcf.gz")
DEFAULT_INTERVALS = Path("/home/leonard.hartmanis/proj/DNA_seq/intervals/panel_targets_all_hg38.interval_list")


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfh) -> None:
    line = f"[{timestamp()}] {msg}"
    print(line)
    logfh.write(line + "\n")
    logfh.flush()


def run_and_log(cmd: list[str], logfh) -> None:
    logfh.write("Command:\n")
    logfh.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
    logfh.flush()
    subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Mutect2 tumor-only calling on one sample.")
    parser.add_argument("--bqsr-dir", default=str(DEFAULT_BQSR_DIR))
    parser.add_argument("--mutect-dir", default=str(DEFAULT_MUTECT_DIR))
    parser.add_argument("--log-dir", default=str(DEFAULT_LOG_DIR))
    parser.add_argument("--reference", default=str(DEFAULT_REF_FA))
    parser.add_argument("--gnomad", default=str(DEFAULT_GNOMAD))
    parser.add_argument("--pon", default=str(DEFAULT_PON))
    parser.add_argument("--intervals", default=str(DEFAULT_INTERVALS))
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--java-options", default="-Xmx4g")
    parser.add_argument("--disable-pon", action="store_true", help="Run without panel of normals")
    parser.add_argument("--disable-intervals", action="store_true", help="Run without -L intervals")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    bqsr_dir = Path(args.bqsr_dir)
    mutect_dir = Path(args.mutect_dir)
    log_dir = Path(args.log_dir)
    reference = Path(args.reference)
    gnomad = Path(args.gnomad)
    pon = Path(args.pon)
    intervals = Path(args.intervals)

    mutect_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    sample_id = args.sample_id

    input_bam = bqsr_dir / f"{sample_id}.bqsr.bam"
    input_bai = bqsr_dir / f"{sample_id}.bqsr.bam.bai"
    if not input_bai.exists():
        alt_bai = bqsr_dir / f"{sample_id}.bqsr.bai"
        if alt_bai.exists():
            input_bai = alt_bai

    unfiltered_vcf = mutect_dir / f"{sample_id}.unfiltered.vcf.gz"
    unfiltered_vcf_tbi = mutect_dir / f"{sample_id}.unfiltered.vcf.gz.tbi"
    filtered_vcf = mutect_dir / f"{sample_id}.filtered.vcf.gz"
    filtered_vcf_tbi = mutect_dir / f"{sample_id}.filtered.vcf.gz.tbi"
    stats_file = mutect_dir / f"{sample_id}.unfiltered.vcf.gz.stats"

    sample_log = log_dir / f"{sample_id}.mutect2.log"

    required = [bqsr_dir, reference, gnomad, input_bam, input_bai]
    if not args.disable_intervals:
        required.append(intervals)
    if not args.disable_pon:
        required.append(pon)

    for p in required:
        if not p.exists():
            raise FileNotFoundError(f"Missing required file/path: {p}")

    expected_outputs = [unfiltered_vcf, filtered_vcf, stats_file]
    if all(p.exists() for p in expected_outputs) and not args.force:
        with open(sample_log, "a") as logfh:
            log(f"SKIP {sample_id}: outputs already exist", logfh)
        return

    mutect_cmd = [
        "gatk", "--java-options", args.java_options,
        "Mutect2",
        "-R", str(reference),
        "-I", str(input_bam),
        "-tumor", sample_id,
        "--germline-resource", str(gnomad),
        "-O", str(unfiltered_vcf),
    ]

    if not args.disable_pon:
        mutect_cmd += ["--panel-of-normals", str(pon)]

    if not args.disable_intervals:
        mutect_cmd += ["-L", str(intervals)]

    filter_cmd = [
        "gatk", "--java-options", args.java_options,
        "FilterMutectCalls",
        "-R", str(reference),
        "-V", str(unfiltered_vcf),
        "--stats", str(stats_file),
        "-O", str(filtered_vcf),
    ]

    cmds = [
        ("Mutect2", mutect_cmd),
        ("FilterMutectCalls", filter_cmd),
    ]

    if args.dry_run:
        for label, cmd in cmds:
            print(f"\n# {label}")
            print(" ".join(shlex.quote(x) for x in cmd))
        return

    start = time.time()

    with open(sample_log, "w") as logfh:
        logfh.write(f"Start time: {datetime.now()}\n")
        logfh.write(f"Sample: {sample_id}\n")
        logfh.write(f"Input BAM: {input_bam}\n")
        logfh.write(f"Reference: {reference}\n")
        logfh.write(f"gnomAD: {gnomad}\n")
        logfh.write(f"PoN: {'DISABLED' if args.disable_pon else pon}\n")
        logfh.write(f"Intervals: {'DISABLED' if args.disable_intervals else intervals}\n")
        logfh.write(f"Java options: {args.java_options}\n\n")
        logfh.flush()

        try:
            for label, cmd in cmds:
                log(f"RUN {label}", logfh)
                run_and_log(cmd, logfh)

            # sanity checks
            for p in [unfiltered_vcf, filtered_vcf, stats_file]:
                if not p.exists():
                    raise RuntimeError(f"Expected output missing: {p}")

            # tabix index can be .tbi or sometimes CSI depending on tool behavior, but usually .tbi here
            if not (unfiltered_vcf_tbi.exists() or Path(str(unfiltered_vcf) + ".csi").exists()):
                raise RuntimeError(f"Missing index for {unfiltered_vcf}")
            if not (filtered_vcf_tbi.exists() or Path(str(filtered_vcf) + ".csi").exists()):
                raise RuntimeError(f"Missing index for {filtered_vcf}")

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
