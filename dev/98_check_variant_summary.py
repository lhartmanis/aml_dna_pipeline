#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

summary_path = Path("/home/leonard.hartmanis/proj/DNA_seq/results/analysis/sample_mutect2_summary_raw.tsv")

df = pd.read_csv(summary_path, sep="\t")

print("\nSamples:", df.shape[0])
print("\nPASS total summary:")
print(df["n_pass_total"].describe())

print("\nTop 20 by PASS count:")
print(df.sort_values("n_pass_total", ascending=False).head(20)[["sample_id", "n_pass_total", "n_pass_snv", "n_pass_indel"]])

print("\nBottom 20 by PASS count:")
print(df.sort_values("n_pass_total", ascending=True).head(20)[["sample_id", "n_pass_total", "n_pass_snv", "n_pass_indel"]])
