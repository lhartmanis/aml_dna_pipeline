#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_QC_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/qc")


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect per-sample BAM QC JSONs into one TSV.")
    parser.add_argument("--qc-dir", default=str(DEFAULT_QC_DIR))
    parser.add_argument("--output", default=None, help="Output TSV path")
    args = parser.parse_args()

    qc_dir = Path(args.qc_dir)
    if not qc_dir.exists():
        raise FileNotFoundError(f"QC dir not found: {qc_dir}")

    output = Path(args.output) if args.output else qc_dir / "bam_qc_summary.tsv"

    rows = []
    for json_file in sorted(qc_dir.glob("*.qc_summary.json")):
        with open(json_file) as fh:
            rows.append(json.load(fh))

    if not rows:
        raise RuntimeError(f"No *.qc_summary.json files found in {qc_dir}")

    df = pd.DataFrame(rows)
    df = df.sort_values("sample_id").reset_index(drop=True)
    df.to_csv(output, sep="\t", index=False)

    print(f"Wrote {len(df)} samples to {output}")


if __name__ == "__main__":
    main()
