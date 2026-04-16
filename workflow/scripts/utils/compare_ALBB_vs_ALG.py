#!/usr/bin/env python3

"""
compare_ALBB_vs_ALG.py

Compare VEP-aware mutation summary metrics between ALBB and ALG cohorts.

What it does:
1. Reads the per-sample summary TSV
2. Infers cohort from sample ID if no explicit cohort/group column exists
3. Produces:
   - summary statistics table
   - Wilcoxon rank-sum (Mann-Whitney U) tests
   - boxplots
   - optional violin plots
4. Focuses on key burden/QC metrics:
   - n_pass_total
   - n_pass_nonsyn
   - n_pass_nonsyn_rare_1e_4
   - burden_total_perMb
   - burden_nonsyn_perMb
   - burden_nonsyn_rare_1e_4_perMb
   - frac_pass
   - vaf_median_pass
   - vaf_mean_pass
   - n_flagged_low_pass
   - n_flagged_low_frac_pass
   - n_flagged_low_rare_nonsyn
   - n_flagged_many_low_vaf

Usage:
python compare_ALBB_vs_ALG.py \
    --input /path/to/sample_vep_qc_flags.tsv \
    --outdir /path/to/vep_cohort_compare

Notes:
- The script tries to be robust to column naming differences.
- If your file already contains a cohort/group column, it will use it.
- Otherwise it will infer cohort from sample IDs beginning with ALBB or ALG.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to sample_vep_qc_flags.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--sample-col",
        default=None,
        help="Optional explicit sample column name. If omitted, script tries to infer it."
    )
    parser.add_argument(
        "--cohort-col",
        default=None,
        help="Optional explicit cohort/group column name. If omitted, script tries to infer it."
    )
    parser.add_argument(
        "--make-violins",
        action="store_true",
        help="Also make violin plots in addition to boxplots."
    )
    return parser.parse_args()


def infer_sample_col(df: pd.DataFrame, user_sample_col: Optional[str] = None) -> str:
    if user_sample_col is not None:
        if user_sample_col not in df.columns:
            raise ValueError(f"--sample-col '{user_sample_col}' not found in columns")
        return user_sample_col

    preferred = ["sample", "Sample", "sample_id", "Sample_ID", "sampleID", "id", "ID"]
    for col in preferred:
        if col in df.columns:
            return col

    # fallback: first object column with many values starting with ALBB/ALG
    object_cols = df.select_dtypes(include=["object"]).columns.tolist()
    for col in object_cols:
        vals = df[col].astype(str)
        frac = ((vals.str.startswith("ALBB")) | (vals.str.startswith("ALG"))).mean()
        if frac > 0.5:
            return col

    raise ValueError("Could not infer sample column. Please provide --sample-col.")


def infer_cohort(
    df: pd.DataFrame,
    sample_col: str,
    user_cohort_col: Optional[str] = None
) -> pd.Series:
    if user_cohort_col is not None:
        if user_cohort_col not in df.columns:
            raise ValueError(f"--cohort-col '{user_cohort_col}' not found in columns")
        cohort = df[user_cohort_col].astype(str)
        return cohort

    for col in ["cohort", "Cohort", "group", "Group", "batch_group"]:
        if col in df.columns:
            return df[col].astype(str)

    sample_ids = df[sample_col].astype(str)
    cohort = np.where(
        sample_ids.str.startswith("ALBB"),
        "ALBB",
        np.where(sample_ids.str.startswith("ALG"), "ALG", "UNKNOWN")
    )
    return pd.Series(cohort, index=df.index, name="cohort")


def choose_metric_columns(df: pd.DataFrame) -> List[str]:
    wanted = [
        "n_pass_total",
        "n_pass_nonsyn",
        "n_pass_nonsyn_rare_1e_4",
        "burden_total_perMb",
        "burden_nonsyn_perMb",
        "burden_nonsyn_rare_1e_4_perMb",
        "frac_pass",
        "vaf_median_pass",
        "vaf_mean_pass",
        "n_flagged_low_pass",
        "n_flagged_low_frac_pass",
        "n_flagged_low_rare_nonsyn",
        "n_flagged_many_low_vaf",
        "n_input_total",
        "n_pass_snv",
        "n_pass_indel",
        "n_pass_missense",
        "n_pass_truncating",
        "n_pass_splice",
    ]
    return [c for c in wanted if c in df.columns]


def safe_numeric(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def summarize_by_cohort(df: pd.DataFrame, metrics: List[str], cohort_col: str = "cohort") -> pd.DataFrame:
    rows = []
    for metric in metrics:
        for cohort, sub in df.groupby(cohort_col):
            vals = sub[metric].dropna()
            if vals.empty:
                continue
            rows.append({
                "metric": metric,
                "cohort": cohort,
                "n": vals.shape[0],
                "mean": vals.mean(),
                "median": vals.median(),
                "std": vals.std(ddof=1) if vals.shape[0] > 1 else np.nan,
                "min": vals.min(),
                "q25": vals.quantile(0.25),
                "q75": vals.quantile(0.75),
                "max": vals.max(),
            })
    return pd.DataFrame(rows)


def mannwhitney_by_metric(df: pd.DataFrame, metrics: List[str], cohort_col: str = "cohort") -> pd.DataFrame:
    rows = []
    sub = df[df[cohort_col].isin(["ALBB", "ALG"])].copy()

    for metric in metrics:
        x = sub.loc[sub[cohort_col] == "ALBB", metric].dropna()
        y = sub.loc[sub[cohort_col] == "ALG", metric].dropna()

        if len(x) == 0 or len(y) == 0:
            continue

        try:
            stat, p = mannwhitneyu(x, y, alternative="two-sided")
        except ValueError:
            stat, p = np.nan, np.nan

        rows.append({
            "metric": metric,
            "n_ALBB": len(x),
            "n_ALG": len(y),
            "median_ALBB": x.median(),
            "median_ALG": y.median(),
            "mean_ALBB": x.mean(),
            "mean_ALG": y.mean(),
            "delta_median_ALG_minus_ALBB": y.median() - x.median(),
            "fold_median_ALG_over_ALBB": (y.median() / x.median()) if x.median() not in [0, np.nan] else np.nan,
            "mannwhitney_u": stat,
            "p_value": p,
        })

    res = pd.DataFrame(rows)
    if not res.empty:
        res["p_adj_bh"] = benjamini_hochberg(res["p_value"].values)
    return res.sort_values("p_value", na_position="last")


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = np.empty(n, dtype=float)

    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev

    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def save_basic_counts(df: pd.DataFrame, outdir: Path, cohort_col: str = "cohort") -> None:
    counts = df[cohort_col].value_counts(dropna=False).rename_axis("cohort").reset_index(name="n_samples")
    counts.to_csv(outdir / "cohort_counts.tsv", sep="\t", index=False)


def make_boxplots(df: pd.DataFrame, metrics: List[str], outdir: Path, cohort_col: str = "cohort") -> None:
    plot_dir = outdir / "boxplots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_df = df[df[cohort_col].isin(["ALBB", "ALG"])].copy()

    for metric in metrics:
        vals = plot_df[[cohort_col, metric]].dropna()
        if vals.empty:
            continue

        albb = vals.loc[vals[cohort_col] == "ALBB", metric].values
        alg = vals.loc[vals[cohort_col] == "ALG", metric].values
        if len(albb) == 0 or len(alg) == 0:
            continue

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.boxplot([albb, alg], labels=["ALBB", "ALG"])
        ax.set_title(metric)
        ax.set_ylabel(metric)

        med_albb = np.median(albb)
        med_alg = np.median(alg)
        ax.text(
            0.5, 0.98,
            f"median ALBB={med_albb:.3g} | median ALG={med_alg:.3g}",
            transform=ax.transAxes,
            ha="center", va="top"
        )

        plt.tight_layout()
        fig.savefig(plot_dir / f"{metric}.png", dpi=200)
        plt.close(fig)


def make_violinplots(df: pd.DataFrame, metrics: List[str], outdir: Path, cohort_col: str = "cohort") -> None:
    plot_dir = outdir / "violins"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_df = df[df[cohort_col].isin(["ALBB", "ALG"])].copy()

    for metric in metrics:
        vals = plot_df[[cohort_col, metric]].dropna()
        if vals.empty:
            continue

        albb = vals.loc[vals[cohort_col] == "ALBB", metric].values
        alg = vals.loc[vals[cohort_col] == "ALG", metric].values
        if len(albb) < 2 or len(alg) < 2:
            continue

        fig, ax = plt.subplots(figsize=(5, 5))
        parts = ax.violinplot([albb, alg], positions=[1, 2], showmeans=True, showmedians=True)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["ALBB", "ALG"])
        ax.set_title(metric)
        ax.set_ylabel(metric)

        plt.tight_layout()
        fig.savefig(plot_dir / f"{metric}.png", dpi=200)
        plt.close(fig)


def write_readme(
    outdir: Path,
    input_path: str,
    sample_col: str,
    metrics: List[str],
) -> None:
    text = f"""Comparison of ALBB vs ALG from VEP summary table

Input:
{input_path}

Detected sample column:
{sample_col}

Metrics analyzed:
{chr(10).join(metrics)}

Outputs:
- cohort_counts.tsv
- summary_by_cohort.tsv
- stats_mannwhitney.tsv
- boxplots/
- violins/   (if requested)
"""
    (outdir / "README.txt").write_text(text)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    sample_col = infer_sample_col(df, args.sample_col)
    df["cohort"] = infer_cohort(df, sample_col, args.cohort_col)

    metrics = choose_metric_columns(df)
    if not metrics:
        raise ValueError("No expected metric columns found in input file.")

    df = safe_numeric(df, metrics)

    # Keep only ALBB / ALG for cohort comparison outputs, but preserve UNKNOWN for counts table if present
    save_basic_counts(df, outdir, cohort_col="cohort")

    summary = summarize_by_cohort(df[df["cohort"].isin(["ALBB", "ALG"])], metrics, cohort_col="cohort")
    summary.to_csv(outdir / "summary_by_cohort.tsv", sep="\t", index=False)

    stats = mannwhitney_by_metric(df[df["cohort"].isin(["ALBB", "ALG"])], metrics, cohort_col="cohort")
    stats.to_csv(outdir / "stats_mannwhitney.tsv", sep="\t", index=False)

    make_boxplots(df, metrics, outdir, cohort_col="cohort")
    if args.make_violins:
        make_violinplots(df, metrics, outdir, cohort_col="cohort")

    write_readme(outdir, args.input, sample_col, metrics)

    print("Done.")
    print(f"Input:  {args.input}")
    print(f"Outdir: {outdir}")
    print(f"Sample column: {sample_col}")
    print(f"Metrics analyzed ({len(metrics)}):")
    for m in metrics:
        print(f"  - {m}")


if __name__ == "__main__":
    main()
