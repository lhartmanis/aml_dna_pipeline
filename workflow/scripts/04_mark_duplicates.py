#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path


DEFAULT_SAMPLE_BAM_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/sample_bam")
DEFAULT_QC_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/qc")
DEFAULT_LOG_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/logs")


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfh) -> None:
    line = f"[{timestamp()}] {msg}"
    print(line)
    logfh.write(line + "\n")
    logfh.flush()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Picard MarkDuplicates on merged sample BAMs.")
    parser.add_argument("--sample-bam-dir", default=str(DEFAULT_SAMPLE_BAM_DIR), help="Directory with merged sample BAMs")
    parser.add_argument("--qc-dir", default=str(DEFAULT_QC_DIR), help="Directory for metrics files")
    parser.add_argument("--log-dir", default=str(DEFAULT_LOG_DIR), help="Directory for logs")
    parser.add_argument("--sample-id", help="Run only one sample")
    parser.add_argument("--java-options", default="-Xmx8g", help='Java options for Picard, e.g. "-Xmx8g"')
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions only")
    args = parser.parse_args()

    sample_bam_dir = Path(args.sample_bam_dir)
    qc_dir = Path(args.qc_dir)
    log_dir = Path(args.log_dir)

    if not sample_bam_dir.exists():
        raise FileNotFoundError(f"Sample BAM directory not found: {sample_bam_dir}")

    qc_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if args.sample_id:
        samples = [args.sample_id]
    else:
        samples = sorted(
            p.name.replace(".merged.bam", "")
            for p in sample_bam_dir.glob("*.merged.bam")
        )

    master_log = log_dir / "mark_duplicates.master.log"
    failures_tsv = log_dir / "mark_duplicates.failures.tsv"
    done_tsv = log_dir / "mark_duplicates.done.tsv"

    with open(master_log, "a") as master_fh, \
         open(failures_tsv, "a") as fail_fh, \
         open(done_tsv, "a") as done_fh:

        if failures_tsv.stat().st_size == 0:
            fail_fh.write("timestamp\tsample_id\treason\n")
            fail_fh.flush()
        if done_tsv.stat().st_size == 0:
            done_fh.write("timestamp\tsample_id\toutput_bam\tmetrics_file\telapsed_sec\n")
            done_fh.flush()

        log(f"START mark_duplicates.py on {len(samples)} samples", master_fh)

        for sample_id in samples:
            start = time.time()

            input_bam = sample_bam_dir / f"{sample_id}.merged.bam"
            input_bai = sample_bam_dir / f"{sample_id}.merged.bam.bai"
            output_bam = sample_bam_dir / f"{sample_id}.markdup.bam"
            output_bai = sample_bam_dir / f"{sample_id}.markdup.bam.bai"
            metrics_file = qc_dir / f"{sample_id}.markdup.metrics.txt"
            sample_log = log_dir / f"{sample_id}.markdup.log"

            if not input_bam.exists():
                reason = f"missing_input_bam:{input_bam}"
                log(f"FAIL {sample_id}: {reason}", master_fh)
                fail_fh.write(f"{timestamp()}\t{sample_id}\t{reason}\n")
                fail_fh.flush()
                continue

            if not input_bai.exists():
                reason = f"missing_input_bai:{input_bai}"
                log(f"FAIL {sample_id}: {reason}", master_fh)
                fail_fh.write(f"{timestamp()}\t{sample_id}\t{reason}\n")
                fail_fh.flush()
                continue

            if (output_bam.exists() or output_bai.exists() or metrics_file.exists()) and not args.force:
                log(f"SKIP {sample_id}: output exists", master_fh)
                continue

            picard_cmd = [
                "picard",
                "-Xmx8g" if args.java_options == "-Xmx8g" else args.java_options,
                "MarkDuplicates",
                f"I={input_bam}",
                f"O={output_bam}",
                f"M={metrics_file}",
                "CREATE_INDEX=true",
                "VALIDATION_STRINGENCY=LENIENT",
                "ASSUME_SORT_ORDER=coordinate",
            ]

            # avoid duplicate -Xmx handling if user gives custom opts
            if args.java_options != "-Xmx8g":
                picard_cmd = [
                    "picard",
                    args.java_options,
                    "MarkDuplicates",
                    f"I={input_bam}",
                    f"O={output_bam}",
                    f"M={metrics_file}",
                    "CREATE_INDEX=true",
                    "VALIDATION_STRINGENCY=LENIENT",
                    "ASSUME_SORT_ORDER=coordinate",
                ]

            if args.dry_run:
                log(f"DRY-RUN {sample_id}", master_fh)
                print(" ".join(shlex.quote(x) for x in picard_cmd))
                continue

            try:
                with open(sample_log, "w") as sample_fh:
                    sample_fh.write(f"Start time: {datetime.now()}\n")
                    sample_fh.write(f"Sample: {sample_id}\n")
                    sample_fh.write(f"Input BAM: {input_bam}\n")
                    sample_fh.write(f"Output BAM: {output_bam}\n")
                    sample_fh.write(f"Metrics: {metrics_file}\n\n")
                    sample_fh.write("Command:\n")
                    sample_fh.write(" ".join(shlex.quote(x) for x in picard_cmd) + "\n\n")
                    sample_fh.flush()

                    subprocess.run(picard_cmd, check=True, stderr=sample_fh, stdout=sample_fh)
                    subprocess.run(
                        ["samtools", "index", str(output_bam)],
                        check=True,
                        stderr=sample_fh,
                        stdout=sample_fh,
                    )
                    end = time.time()
                    elapsed = end - start

                    sample_fh.write(f"\nEnd time: {datetime.now()}\n")
                    sample_fh.write(f"Elapsed seconds: {elapsed:.2f}\n")
                    sample_fh.flush()

                    picard_bai = output_bam.with_suffix(".bai")          # ALBB13002.markdup.bai
                    samtools_bai = Path(str(output_bam) + ".bai")        # ALBB13002.markdup.bam.bai

                    if not output_bam.exists():
                        raise RuntimeError(f"Output BAM missing: {output_bam}")
                    if not picard_bai.exists() and not samtools_bai.exists():
                        raise RuntimeError(
                            f"Output BAI missing: checked {picard_bai} and {samtools_bai}"
                        )
                    if not metrics_file.exists():
                        raise RuntimeError(f"Metrics file missing: {metrics_file}")

                log(f"DONE {sample_id}: elapsed={elapsed:.2f}s", master_fh)
                done_fh.write(f"{timestamp()}\t{sample_id}\t{output_bam}\t{metrics_file}\t{elapsed:.2f}\n")
                done_fh.flush()

            except Exception as e:
                log(f"FAIL {sample_id}: {e}", master_fh)
                fail_fh.write(f"{timestamp()}\t{sample_id}\t{str(e)}\n")
                fail_fh.flush()

        log("END mark_duplicates.py", master_fh)


if __name__ == "__main__":
    main()
