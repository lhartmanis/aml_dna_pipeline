#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfh) -> None:
    line = f"[{timestamp()}] {msg}"
    print(line)
    logfh.write(line + "\n")
    logfh.flush()


def run_and_log(cmd: list[str], logfh, env: dict[str, str] | None = None) -> None:
    logfh.write("Command:\n")
    logfh.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
    logfh.flush()
    subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh, env=env)


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate one normalized VCF with Ensembl VEP.")
    parser.add_argument("--input-dir", required=True, help="Directory containing normalized VCFs")
    parser.add_argument("--output-dir", required=True, help="Directory for VEP-annotated VCFs")
    parser.add_argument("--log-dir", required=True, help="Directory for per-sample logs")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--vep-cache", required=True, help="VEP cache root directory")
    parser.add_argument("--vep-fasta", required=True, help="FASTA for offline VEP annotation")
    parser.add_argument("--assembly", default="GRCh38", help="Genome assembly")
    parser.add_argument("--vep-exe", default="vep", help="VEP executable")
    parser.add_argument("--force", action="store_true", help="Overwrite outputs if they exist")
    parser.add_argument("--dry-run", action="store_true", help="Print commands only")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    log_dir = Path(args.log_dir)
    vep_cache = Path(args.vep_cache)
    vep_fasta = Path(args.vep_fasta)

    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    sample_id = args.sample_id

    input_vcf = input_dir / f"{sample_id}.filtered.norm.vcf.gz"
    input_tbi = input_dir / f"{sample_id}.filtered.norm.vcf.gz.tbi"

    output_vcf = output_dir / f"{sample_id}.filtered.vep.vcf.gz"
    output_tbi = output_dir / f"{sample_id}.filtered.vep.vcf.gz.tbi"

    sample_log = log_dir / f"{sample_id}.vep.log"

    if not input_vcf.exists():
        raise FileNotFoundError(f"Missing input VCF: {input_vcf}")
    if not input_tbi.exists():
        raise FileNotFoundError(f"Missing input VCF index: {input_tbi}")
    if not vep_cache.exists():
        raise FileNotFoundError(f"Missing VEP cache directory: {vep_cache}")
    if not vep_fasta.exists():
        raise FileNotFoundError(f"Missing VEP FASTA: {vep_fasta}")

    if output_vcf.exists() and output_tbi.exists() and not args.force:
        with open(sample_log, "a") as logfh:
            log(f"SKIP {sample_id}: outputs already exist", logfh)
        return

    vep_cmd = [
        args.vep_exe,
        "--cache",
        "--offline",
        "--assembly", args.assembly,
        "--dir_cache", str(vep_cache),
        "--fasta", str(vep_fasta),
        "--input_file", str(input_vcf),
        "--output_file", str(output_vcf),
        "--vcf",
        "--compress_output", "bgzip",
        "--force_overwrite",
        "--symbol",
        "--hgvs",
        "--terms", "SO",
        "--pick",
        "--pick_order", "canonical,appris,biotype,rank,length",
        "--everything",
    ]

    tabix_cmd = [
        "tabix",
        "-f",
        "-p", "vcf",
        str(output_vcf),
    ]

    if args.dry_run:
        print("# VEP")
        print(" ".join(shlex.quote(x) for x in vep_cmd))
        print("\n# tabix")
        print(" ".join(shlex.quote(x) for x in tabix_cmd))
        return

    start = time.time()

    with open(sample_log, "w") as logfh:
        logfh.write(f"Start time: {datetime.now()}\n")
        logfh.write(f"Sample: {sample_id}\n")
        logfh.write(f"Input VCF: {input_vcf}\n")
        logfh.write(f"Output VCF: {output_vcf}\n")
        logfh.write(f"VEP cache: {vep_cache}\n")
        logfh.write(f"VEP FASTA: {vep_fasta}\n")
        logfh.write(f"Assembly: {args.assembly}\n\n")
        logfh.flush()

        try:
            env = None
            log(f"RUN VEP", logfh)
            run_and_log(vep_cmd, logfh, env=env)

            log(f"RUN tabix", logfh)
            run_and_log(tabix_cmd, logfh)

            if not output_vcf.exists():
                raise RuntimeError(f"Expected output missing: {output_vcf}")
            if not output_tbi.exists():
                raise RuntimeError(f"Expected output index missing: {output_tbi}")

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
