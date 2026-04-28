#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path

import pandas as pd


DEFAULT_MANIFEST = Path("results/manifest/fastq_manifest.tsv")
DEFAULT_UNIT_BAM_DIR = Path("results/bam")
DEFAULT_SAMPLE_BAM_DIR = Path("results/sample_bam")
DEFAULT_LOG_DIR = Path("logs/merge")


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfh) -> None:
    line = f"[{timestamp()}] {msg}"
    print(line)
    logfh.write(line + "\n")
    logfh.flush()


def safe_unlink(path: Path) -> None:
    if path.exists() or path.is_symlink():
        path.unlink()


def ensure_symlink(src: Path, dst: Path) -> None:
    safe_unlink(dst)
    os.symlink(src.resolve(), dst)


def samtools_index(bam: Path, logfh) -> None:
    cmd = ["samtools", "index", str(bam)]
    logfh.write("Index command:\n")
    logfh.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
    logfh.flush()
    subprocess.run(cmd, check=True, stderr=logfh)


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge unit BAMs per sample using the manifest.")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="Path to fastq_manifest.tsv")
    parser.add_argument("--unit-bam-dir", default=str(DEFAULT_UNIT_BAM_DIR), help="Directory with unit-level BAMs")
    parser.add_argument("--sample-bam-dir", default=str(DEFAULT_SAMPLE_BAM_DIR), help="Directory for merged sample BAMs")
    parser.add_argument("--log-dir", default=str(DEFAULT_LOG_DIR), help="Directory for logs")
    parser.add_argument("--sample-id", help="Run only one sample")
    parser.add_argument("--samtools-threads", type=int, default=4, help="Threads for samtools merge")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions only")
    args = parser.parse_args()

    manifest = Path(args.manifest)
    unit_bam_dir = Path(args.unit_bam_dir)
    sample_bam_dir = Path(args.sample_bam_dir)
    log_dir = Path(args.log_dir)

    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")
    if not unit_bam_dir.exists():
        raise FileNotFoundError(f"Unit BAM directory not found: {unit_bam_dir}")

    sample_bam_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(manifest, sep="\t")
    if args.sample_id:
        df = df[df["sample_id"] == args.sample_id].copy()
        if df.empty:
            raise ValueError(f"No rows found for sample_id={args.sample_id}")

    samples = sorted(df["sample_id"].unique())

    master_log = log_dir / "merge_sample_bams.master.log"
    failures_tsv = log_dir / "merge_sample_bams.failures.tsv"
    done_tsv = log_dir / "merge_sample_bams.done.tsv"

    with open(master_log, "a") as master_fh, \
         open(failures_tsv, "a") as fail_fh, \
         open(done_tsv, "a") as done_fh:

        if failures_tsv.stat().st_size == 0:
            fail_fh.write("timestamp\tsample_id\treason\n")
            fail_fh.flush()
        if done_tsv.stat().st_size == 0:
            done_fh.write("timestamp\tsample_id\toutput_bam\tmode\telapsed_sec\n")
            done_fh.flush()

        log(f"START merge_sample_bams.py on {len(samples)} samples", master_fh)

        for sample_id in samples:
            start = time.time()
            sample_rows = df[df["sample_id"] == sample_id].copy()

            expected_rgids = sample_rows["rgid"].tolist()
            unit_bams = [unit_bam_dir / f"{rgid}.sorted.bam" for rgid in expected_rgids]
            unit_bais = [unit_bam_dir / f"{rgid}.sorted.bam.bai" for rgid in expected_rgids]

            missing = [str(p) for p in unit_bams if not p.exists()]
            missing_idx = [str(p) for p in unit_bais if not p.exists()]

            if missing or missing_idx:
                reason_parts = []
                if missing:
                    reason_parts.append(f"missing_bams={len(missing)}")
                if missing_idx:
                    reason_parts.append(f"missing_bais={len(missing_idx)}")
                reason = ";".join(reason_parts)
                log(f"FAIL {sample_id}: {reason}", master_fh)
                fail_fh.write(f"{timestamp()}\t{sample_id}\t{reason}\n")
                fail_fh.flush()
                continue

            merged_bam = sample_bam_dir / f"{sample_id}.merged.bam"
            merged_bai = sample_bam_dir / f"{sample_id}.merged.bam.bai"
            sample_log = log_dir / f"{sample_id}.merge.log"

            if (merged_bam.exists() or merged_bai.exists()) and not args.force:
                log(f"SKIP {sample_id}: output exists", master_fh)
                continue

            if args.dry_run:
                mode = "symlink" if len(unit_bams) == 1 else "merge"
                log(f"DRY-RUN {sample_id}: mode={mode}, n_units={len(unit_bams)}", master_fh)
                continue

            try:
                with open(sample_log, "w") as sample_fh:
                    sample_fh.write(f"Start time: {datetime.now()}\n")
                    sample_fh.write(f"Sample: {sample_id}\n")
                    sample_fh.write(f"Expected units: {len(unit_bams)}\n")
                    sample_fh.write("Input BAMs:\n")
                    for bam in unit_bams:
                        sample_fh.write(f"  {bam}\n")
                    sample_fh.write("\n")
                    sample_fh.flush()

                    if len(unit_bams) == 1:
                        src_bam = unit_bams[0]
                        src_bai = unit_bais[0]

                        sample_fh.write("Mode: symlink\n\n")
                        sample_fh.flush()

                        safe_unlink(merged_bam)
                        safe_unlink(merged_bai)

                        ensure_symlink(src_bam, merged_bam)
                        ensure_symlink(src_bai, merged_bai)

                        mode = "symlink"

                    else:
                        sample_fh.write("Mode: merge\n\n")
                        merge_cmd = [
                            "samtools", "merge",
                            "-@", str(args.samtools_threads),
                            "-o", str(merged_bam),
                            *[str(x) for x in unit_bams]
                        ]

                        sample_fh.write("Merge command:\n")
                        sample_fh.write(" ".join(shlex.quote(x) for x in merge_cmd) + "\n\n")
                        sample_fh.flush()

                        safe_unlink(merged_bam)
                        safe_unlink(merged_bai)

                        subprocess.run(merge_cmd, check=True, stderr=sample_fh)

                        samtools_index(merged_bam, sample_fh)

                        mode = "merge"

                    end = time.time()
                    elapsed = end - start

                    sample_fh.write(f"End time: {datetime.now()}\n")
                    sample_fh.write(f"Elapsed seconds: {elapsed:.2f}\n")
                    sample_fh.flush()

                if not merged_bam.exists():
                    raise RuntimeError(f"Output BAM missing after processing: {merged_bam}")
                if not merged_bai.exists():
                    raise RuntimeError(f"Output BAI missing after processing: {merged_bai}")

                log(f"DONE {sample_id}: mode={mode}, n_units={len(unit_bams)}, elapsed={elapsed:.2f}s", master_fh)
                done_fh.write(f"{timestamp()}\t{sample_id}\t{merged_bam}\t{mode}\t{elapsed:.2f}\n")
                done_fh.flush()

            except Exception as e:
                log(f"FAIL {sample_id}: {e}", master_fh)
                fail_fh.write(f"{timestamp()}\t{sample_id}\t{str(e)}\n")
                fail_fh.flush()

        log("END merge_sample_bams.py", master_fh)


if __name__ == "__main__":
    main()
