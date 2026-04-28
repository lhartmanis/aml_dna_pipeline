#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_QC_DIR = Path("results/qc")
DEFAULT_OUTPUT = Path("results/qc/bam_qc_summary.tsv")


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect per-sample BAM QC JSONs into one TSV.")
    parser.add_argument("--qc-dir", default=str(DEFAULT_QC_DIR), help="Directory containing *.qc_summary.json files")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT), help="Output TSV path")
    args = parser.parse_args()

    qc_dir = Path(args.qc_dir)
    output = Path(args.output)

    if not qc_dir.exists():
        raise FileNotFoundError(f"QC dir not found: {qc_dir}")

    rows = []
    for json_file in sorted(qc_dir.glob("*.qc_summary.json")):
        with open(json_file) as fh:
            rows.append(json.load(fh))

    if not rows:
        raise RuntimeError(f"No *.qc_summary.json files found in {qc_dir}")

    df = pd.DataFrame(rows)
    if "sample_id" in df.columns:
        df = df.sort_values("sample_id").reset_index(drop=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)

    print(f"Wrote {len(df)} samples to {output}")


if __name__ == "__main__":
    main()
