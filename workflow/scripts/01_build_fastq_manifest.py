#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path
import pandas as pd

ROOT = Path("/home/leonard.hartmanis/proj/DNA_seq/raw_fastq")
OUTDIR = Path("/home/leonard.hartmanis/proj/DNA_seq/metadata")
MANIFEST_OUT = OUTDIR / "fastq_manifest.tsv"
PROBLEMS_OUT = OUTDIR / "fastq_manifest_problems.tsv"
SUMMARY_OUT = OUTDIR / "fastq_manifest_summary.txt"

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")


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

    # e.g. ..._R1_001 or ..._R2_001
    m = re.search(r"(?:^|[_\.])R([12])(?:[_\.]\d+)$", stem)
    if m:
        return f"R{m.group(1)}"

    # e.g. ..._R1 or ..._R2
    m = re.search(r"(?:^|[_\.])R([12])$", stem)
    if m:
        return f"R{m.group(1)}"

    # e.g. ..._1 or ..._2
    m = re.search(r"(?:[_\.])([12])$", stem)
    if m:
        return f"R{m.group(1)}"

    return None


def strip_read_token(filename: str) -> str:
    """
    Remove the terminal read token only.
    """
    stem = remove_fastq_suffix(filename)

    # Remove ..._R1_001 / ..._R2_001
    stem = re.sub(r"(?:[_\.])R[12](?:[_\.]\d+)$", "", stem)

    # Remove ..._R1 / ..._R2
    stem = re.sub(r"(?:[_\.])R[12]$", "", stem)

    # Remove ..._1 / ..._2
    stem = re.sub(r"(?:[_\.])[12]$", "", stem)

    # Clean trailing separators
    stem = re.sub(r"[_\.]+$", "", stem)
    return stem


def infer_lane(filename: str) -> str:
    """
    Tries several naming conventions.
    """
    stem = remove_fastq_suffix(filename)

    # e.g. ..._L003_...
    m = re.search(r"_L(\d{3})_", stem)
    if m:
        return f"L{m.group(1)}"

    # e.g. lane3
    m = re.search(r"lane(\d+)", stem, flags=re.IGNORECASE)
    if m:
        return f"L{m.group(1)}"

    # e.g. 7_141222_...
    m = re.match(r"^(\d+)_", stem)
    if m:
        return f"L{m.group(1)}"

    return "NA"


def infer_flowcell(filename: str) -> str:
    """
    Capture common Illumina flowcell-like tokens.
    """
    stem = remove_fastq_suffix(filename)

    # e.g. BC6KGRANXX, AC6KH3ANXX, BC6FA6ANXX
    m = re.search(r"(?:^|_)([A-Z0-9]{8,}ANXX)(?:_|$)", stem)
    if m:
        return m.group(1)

    # fallback: any long alnum token ending in ANXX
    m = re.search(r"(?:^|_)([A-Z0-9]+ANXX)(?:_|$)", stem)
    if m:
        return m.group(1)

    return "NA"


def infer_library(sample_id: str) -> str:
    return sample_id


def infer_platform_unit(flowcell: str, lane: str, unit_base: str) -> str:
    if flowcell != "NA" and lane != "NA":
        return f"{flowcell}.{lane}"
    return unit_base


def main() -> None:
    sample_dirs = sorted([p for p in ROOT.iterdir() if p.is_dir()])

    rows: list[dict] = []
    problems: list[dict] = []

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        fastqs = sorted(
            [
                p for p in sample_dir.iterdir()
                if p.is_file() and p.name.endswith(FASTQ_SUFFIXES)
            ]
        )

        if not fastqs:
            problems.append(
                {
                    "sample_id": sample_id,
                    "item": str(sample_dir),
                    "problem": "No FASTQ files found"
                }
            )
            continue

        grouped: dict[str, dict[str, Path]] = {}

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

            if read in grouped[unit_base]:
                problems.append(
                    {
                        "sample_id": sample_id,
                        "item": fq.name,
                        "problem": f"Duplicate {read} for unit {unit_base}"
                    }
                )
            else:
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

    problems_df = pd.DataFrame(problems)
    if not problems_df.empty:
        problems_df = problems_df.sort_values(["sample_id", "item"]).reset_index(drop=True)

    manifest.to_csv(MANIFEST_OUT, sep="\t", index=False)
    problems_df.to_csv(PROBLEMS_OUT, sep="\t", index=False)

    n_sample_dirs = len(sample_dirs)
    n_samples_in_manifest = manifest["sample_id"].nunique() if not manifest.empty else 0
    n_units = len(manifest)
    samples_with_multiple_units = int((manifest.groupby("sample_id").size() > 1).sum()) if not manifest.empty else 0
    n_problem_rows = len(problems_df)

    with open(SUMMARY_OUT, "w") as fh:
        fh.write(f"n_sample_dirs\t{n_sample_dirs}\n")
        fh.write(f"n_samples_in_manifest\t{n_samples_in_manifest}\n")
        fh.write(f"n_units\t{n_units}\n")
        fh.write(f"samples_with_multiple_units\t{samples_with_multiple_units}\n")
        fh.write(f"n_problem_rows\t{n_problem_rows}\n")

    print(f"Manifest written: {MANIFEST_OUT}")
    print(f"Problems written: {PROBLEMS_OUT}")
    print(f"Summary written:  {SUMMARY_OUT}")
    print()
    print(f"n_sample_dirs: {n_sample_dirs}")
    print(f"n_samples_in_manifest: {n_samples_in_manifest}")
    print(f"n_units: {n_units}")
    print(f"samples_with_multiple_units: {samples_with_multiple_units}")
    print(f"n_problem_rows: {n_problem_rows}")


if __name__ == "__main__":
    main()

