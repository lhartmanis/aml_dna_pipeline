#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path

import pandas as pd


MANIFEST = Path("/home/leonard.hartmanis/proj/DNA_seq/metadata/fastq_manifest.tsv")
REF = Path("/home/leonard.hartmanis/proj/DNA_seq/ref/GRCh38.p14.genome.fa")
OUTDIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/unit_bam")
LOGDIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/logs")


def build_rg_string(row: pd.Series) -> str:
    fields = {
        "ID": row["rgid"],
        "SM": row["rgsm"],
        "LB": row["rglb"],
        "PL": row["rgpl"],
        "PU": row["rgpu"],
    }
    return "@RG\\t" + "\\t".join(f"{k}:{v}" for k, v in fields.items())


def select_row(df: pd.DataFrame, sample_id: str | None, unit_base: str | None, rgid: str | None) -> pd.Series:
    if rgid:
        sub = df[df["rgid"] == rgid].copy()
        if len(sub) == 0:
            raise ValueError(f"No row found for rgid={rgid}")
        if len(sub) > 1:
            raise ValueError(f"Multiple rows found for rgid={rgid}, expected exactly one")
        return sub.iloc[0]

    if not sample_id:
        raise ValueError("Provide either --rgid or --sample-id")

    sub = df[df["sample_id"] == sample_id].copy()
    if len(sub) == 0:
        raise ValueError(f"No rows found for sample_id={sample_id}")

    if unit_base:
        sub = sub[sub["unit_base"] == unit_base].copy()
        if len(sub) == 0:
            raise ValueError(f"No row found for sample_id={sample_id} and unit_base={unit_base}")
        if len(sub) > 1:
            raise ValueError(f"Multiple rows found for sample_id={sample_id} and unit_base={unit_base}")
        return sub.iloc[0]

    if len(sub) > 1:
        choices = "\n".join(sub["unit_base"].tolist())
        raise ValueError(
            f"Sample {sample_id} has multiple units. Re-run with --unit-base.\nChoices:\n{choices}"
        )

    return sub.iloc[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="Align one FASTQ pair from the manifest with bwa-mem2 and samtools sort.")
    parser.add_argument("--sample-id", help="Sample ID from manifest")
    parser.add_argument("--unit-base", help="Unit base from manifest for multi-unit samples")
    parser.add_argument("--rgid", help="RGID from manifest; alternative to sample/unit selection")
    parser.add_argument("--manifest", default=str(MANIFEST), help="Path to fastq_manifest.tsv")
    parser.add_argument("--reference", default=str(REF), help="Path to reference FASTA")
    parser.add_argument("--outdir", default=str(OUTDIR), help="Output directory for unit BAMs")
    parser.add_argument("--logdir", default=str(LOGDIR), help="Output directory for logs")
    parser.add_argument("--bwa-threads", type=int, default=8, help="Threads for bwa-mem2 mem")
    parser.add_argument("--sort-threads", type=int, default=4, help="Threads for samtools sort")
    parser.add_argument("--sort-memory", default="1G", help="Memory per thread for samtools sort, e.g. 1G")
    parser.add_argument("--force", action="store_true", help="Overwrite existing BAM and BAI")
    parser.add_argument("--dry-run", action="store_true", help="Print commands and selected row without executing")
    args = parser.parse_args()

    start_time = time.time()
    start_dt = datetime.now()

    manifest = Path(args.manifest)
    reference = Path(args.reference)
    outdir = Path(args.outdir)
    logdir = Path(args.logdir)

    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")

    outdir.mkdir(parents=True, exist_ok=True)
    logdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(manifest, sep="\t")
    row = select_row(df, args.sample_id, args.unit_base, args.rgid)

    fastq_r1 = Path(row["fastq_r1"])
    fastq_r2 = Path(row["fastq_r2"])
    if not fastq_r1.exists():
        raise FileNotFoundError(f"R1 not found: {fastq_r1}")
    if not fastq_r2.exists():
        raise FileNotFoundError(f"R2 not found: {fastq_r2}")

    rgid = row["rgid"]
    bam_out = outdir / f"{rgid}.sorted.bam"
    bai_out = outdir / f"{rgid}.sorted.bam.bai"
    log_out = logdir / f"{rgid}.align.log"

    if bam_out.exists() and not args.force:
        raise FileExistsError(f"Output BAM already exists: {bam_out} (use --force to overwrite)")

    rg_string = build_rg_string(row)

    bwa_cmd = [
        "bwa-mem2", "mem",
        "-t", str(args.bwa_threads),
        "-R", rg_string,
        str(reference),
        str(fastq_r1),
        str(fastq_r2),
    ]

    sort_cmd = [
        "samtools", "sort",
        "-@", str(args.sort_threads),
        "-m", args.sort_memory,
        "-o", str(bam_out),
        "-",
    ]

    index_cmd = ["samtools", "index", str(bam_out)]

    if args.dry_run:
        print("Selected row:")
        print(row.to_string())
        print("\nBWA command:")
        print(" ".join(shlex.quote(x) for x in bwa_cmd))
        print("\nSort command:")
        print(" ".join(shlex.quote(x) for x in sort_cmd))
        print("\nIndex command:")
        print(" ".join(shlex.quote(x) for x in index_cmd))
        print(f"\nLog file: {log_out}")
        return

    with open(log_out, "w") as logfh:
        logfh.write(f"Start time: {start_dt}\n\n")
        logfh.write("BWA command:\n")
        logfh.write(" ".join(shlex.quote(x) for x in bwa_cmd) + "\n\n")
        logfh.write("Sort command:\n")
        logfh.write(" ".join(shlex.quote(x) for x in sort_cmd) + "\n\n")

        bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=logfh, text=False)
        sort_proc = subprocess.Popen(sort_cmd, stdin=bwa_proc.stdout, stderr=logfh, text=False)

        assert bwa_proc.stdout is not None
        bwa_proc.stdout.close()

        sort_return = sort_proc.wait()
        bwa_return = bwa_proc.wait()

        if bwa_return != 0:
            raise RuntimeError(f"bwa-mem2 failed with exit code {bwa_return}. See log: {log_out}")
        if sort_return != 0:
            raise RuntimeError(f"samtools sort failed with exit code {sort_return}. See log: {log_out}")

        logfh.write("\nIndex command:\n")
        logfh.write(" ".join(shlex.quote(x) for x in index_cmd) + "\n")

        subprocess.run(index_cmd, check=True, stderr=logfh)

        end_time = time.time()
        end_dt = datetime.now()
        elapsed_sec = end_time - start_time
        elapsed_min = elapsed_sec / 60

        logfh.write(f"\nEnd time: {end_dt}\n")
        logfh.write(f"Elapsed seconds: {elapsed_sec:.2f}\n")
        logfh.write(f"Elapsed minutes: {elapsed_min:.2f}\n")

    if not bam_out.exists():
        raise RuntimeError(f"Expected BAM not found after alignment: {bam_out}")
    if not bai_out.exists():
        raise RuntimeError(f"Expected BAM index not found after indexing: {bai_out}")

    print(f"Done: {bam_out}")
    print(f"Index: {bai_out}")
    print(f"Log:   {log_out}")
    print("\n=== Alignment timing ===")
    print(f"Start:   {start_dt}")
    print(f"End:     {end_dt}")
    print(f"Elapsed: {elapsed_sec:.2f} sec ({elapsed_min:.2f} min)")



if __name__ == "__main__":
    main()
