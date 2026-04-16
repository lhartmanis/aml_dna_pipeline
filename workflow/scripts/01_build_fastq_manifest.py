#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
import pandas as pd

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build FASTQ manifest for AML DNA panel pipeline."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Root directory containing one subdirectory per sample with FASTQ files."
    )
    parser.add_argument(
        "--manifest-out",
        required=True,
        help="Output TSV path for the FASTQ manifest."
    )
    parser.add_argument(
        "--problems-out",
        required=True,
        help="Output TSV path for FASTQ pairing / parsing problems."
    )
    parser.add_argument(
        "--summary-out",
        required=True,
        help="Output TXT path for manifest summary statistics."
    )
    return parser.parse_args()


def remove_fastq_suffix(filename: str) -> str:
    for suf in FASTQ_SUFFIXES:
        if filename.endswith(suf):
            return filename[: -len(suf)]
    return filename


def detect_read(filename: str) -> str | None:
    """
    Detect read direction from the END of the filename stem only.

    Supported:
      sample_R1.fastq.gz
      sample_R2.fastq.gz
      sample_R1_001.fastq.gz
      sample_R2_001.fastq.gz
      sample_1.fastq.gz
      sample_2.fastq.gz
    """
    stem = remove_fastq_suffix(filename)

    m = re.search(r"(?:^|[_\.])R([12])(?:[_\.]\d+)$", stem)
    if m:
        return f"R{m.group(1)}"

    m = re.search(r"(?:^|[_\.])R([12])$", stem)
    if m:
        return f"R{m.group(1)}"

    m = re.search(r"(?:[_\.])([12])$", stem)
    if m:
        return f"R{m.group(1)}"

    return None


def strip_read_token(filename: str) -> str:
    stem = remove_fastq_suffix(filename)
    stem = re.sub(r"(?:[_\.])R[12](?:[_\.]\d+)$", "", stem)
    stem = re.sub(r"(?:[_\.])R[12]$", "", stem)
    stem = re.sub(r"(?:[_\.])[12]$", "", stem)
    stem = re.sub(r"[_\.]+$", "", stem)
    return stem


def infer_lane(filename: str) -> str:
    stem = remove_fastq_suffix(filename)

    m = re.search(r"(?:^|[_\.])L(\d{3})(?:[_\.]|$)", stem)
    if m:
        return f"L{m.group(1)}"

    m = re.search(r"(?:^|[_\.])lane(\d+)(?:[_\.]|$)", stem, flags=re.IGNORECASE)
    if m:
        return f"L{int(m.group(1)):03d}"

    return "NA"


def infer_flowcell(filename: str) -> str:
    stem = remove_fastq_suffix(filename)

    # Example: 7_141222_BC6KGRANXX_P1781_1336_1.fastq.gz
    m = re.search(r"\d+_(\d{6})_([A-Z0-9]+XX)(?:_|$)", stem)
    if m:
        return m.group(2)

    m = re.search(r"([A-Z0-9]+XX)", stem)
    if m:
        return m.group(1)

    return "NA"


def infer_library(sample_id: str) -> str:
    return sample_id


def infer_platform_unit(flowcell: str, lane: str, unit_base: str) -> str:
    return f"{flowcell}.{lane}" if flowcell != "NA" else f"NA.{lane}"


def main():
    args = parse_args()

    root = Path(args.input_dir)
    manifest_out = Path(args.manifest_out)
    problems_out = Path(args.problems_out)
    summary_out = Path(args.summary_out)

    if not root.exists():
        raise FileNotFoundError(f"Input directory does not exist: {root}")

    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    problems_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted([p for p in root.iterdir() if p.is_dir()])

    rows = []
    problems = []

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        fastqs = sorted(
            [p for p in sample_dir.iterdir() if p.is_file() and p.name.endswith(FASTQ_SUFFIXES)]
        )

        grouped = {}
        for fq in fastqs:
            read = detect_read(fq.name)
            if read is None:
                problems.append(
                    {
                        "sample_id": sample_id,
                        "item": fq.name,
                        "problem": "Could not determine read direction"
                    }
                )
                continue

            unit_base = strip_read_token(fq.name)
            grouped.setdefault(unit_base, {})
            grouped[unit_base][read] = fq

        for unit_base in sorted(grouped):
            pair = grouped[unit_base]
            r1 = pair.get("R1")
            r2 = pair.get("R2")

            if r1 is None or r2 is None:
                missing = []
                if r1 is None:
                    missing.append("R1")
                if r2 is None:
                    missing.append("R2")

                problems.append(
                    {
                        "sample_id": sample_id,
                        "item": unit_base,
                        "problem": f"Missing mate: {','.join(missing)}"
                    }
                )
                continue

            lane = infer_lane(r1.name)
            flowcell = infer_flowcell(r1.name)
            library = infer_library(sample_id)

            rgsm = sample_id
            rglb = library
            rgpl = "ILLUMINA"
            rgpu = infer_platform_unit(flowcell, lane, unit_base)
            rgid = f"{sample_id}.{rgpu}"

            rows.append(
                {
                    "sample_id": sample_id,
                    "unit_base": unit_base,
                    "fastq_r1": str(r1),
                    "fastq_r2": str(r2),
                    "lane": lane,
                    "flowcell": flowcell,
                    "library": library,
                    "platform": rgpl,
                    "rgid": rgid,
                    "rgsm": rgsm,
                    "rglb": rglb,
                    "rgpl": rgpl,
                    "rgpu": rgpu,
                }
            )

    manifest = pd.DataFrame(rows)
    if not manifest.empty:
        manifest = manifest.sort_values(["sample_id", "unit_base"]).reset_index(drop=True)

    problems_df = pd.DataFrame(
        problems,
        columns=["sample_id", "item", "problem"]
    )
    if not problems_df.empty:
        problems_df = problems_df.sort_values(["sample_id", "item"]).reset_index(drop=True)

    manifest.to_csv(manifest_out, sep="\t", index=False)
    problems_df.to_csv(problems_out, sep="\t", index=False)

    n_sample_dirs = len(sample_dirs)
    n_samples_in_manifest = manifest["sample_id"].nunique() if not manifest.empty else 0
    n_units = len(manifest)
    samples_with_multiple_units = int((manifest.groupby("sample_id").size() > 1).sum()) if not manifest.empty else 0
    n_problem_rows = len(problems_df)

    with open(summary_out, "w") as fh:
        fh.write(f"n_sample_dirs\t{n_sample_dirs}\n")
        fh.write(f"n_samples_in_manifest\t{n_samples_in_manifest}\n")
        fh.write(f"n_units\t{n_units}\n")
        fh.write(f"samples_with_multiple_units\t{samples_with_multiple_units}\n")
        fh.write(f"n_problem_rows\t{n_problem_rows}\n")

    print(f"Manifest written: {manifest_out}")
    print(f"Problems written: {problems_out}")
    print(f"Summary written:  {summary_out}")


if __name__ == "__main__":
    main()
