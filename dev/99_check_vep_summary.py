#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq")
ANALYSIS_DIR = PROJECT_DIR / "results" / "analysis"
OUT_DIR = ANALYSIS_DIR / "qc_vep"

SUMMARY_FILE = ANALYSIS_DIR / "sample_mutect2_summary_annotated.tsv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

QC_TABLE = OUT_DIR / "sample_vep_qc_flags.tsv"
QC_REPORT = OUT_DIR / "vep_qc_report.txt"


KEY_METRICS = [
    "n_total",
    "n_pass_total",
    "n_pass_nonsyn",
    "n_pass_trunc",
    "n_pass_splice",
    "n_pass_nonsyn_rare_1e_3",
    "n_pass_nonsyn_rare_1e_4",
    "burden_total_perMb",
    "burden_nonsyn_perMb",
    "burden_nonsyn_rare_1e_3_perMb",
    "burden_nonsyn_rare_1e_4_perMb",
    "frac_pass",
    "vaf_median_pass",
    "vaf_max_pass",
    "frac_vaf_ge_0_25",
    "frac_vaf_ge_0_40",
    "n_low_vaf_lt_0_05",
    "dp_median_pass",
    "dp_iqr_pass",
    "ad_alt_median_pass",
    "frac_pass_with_existing_variation",
    "frac_pass_with_rsid",
    "frac_pass_with_max_af",
    "frac_pass_with_max_af_gt_0_01",
]


def robust_z(x: pd.Series) -> pd.Series:
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if np.isnan(mad) or mad == 0:
        return pd.Series(np.nan, index=x.index)
    return 0.6745 * (x - med) / mad


def summarize_series(s: pd.Series) -> dict[str, float]:
    vals = s.dropna()
    if len(vals) == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "median": np.nan,
            "sd": np.nan,
            "min": np.nan,
            "p05": np.nan,
            "p25": np.nan,
            "p75": np.nan,
            "p95": np.nan,
            "max": np.nan,
        }
    return {
        "n": int(vals.shape[0]),
        "mean": float(vals.mean()),
        "median": float(vals.median()),
        "sd": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
        "min": float(vals.min()),
        "p05": float(vals.quantile(0.05)),
        "p25": float(vals.quantile(0.25)),
        "p75": float(vals.quantile(0.75)),
        "p95": float(vals.quantile(0.95)),
        "max": float(vals.max()),
    }


def add_flag(df: pd.DataFrame, col: str, mask: pd.Series, reason: str) -> None:
    df[reason] = mask.astype(int)


def make_hist(df: pd.DataFrame, col: str, out_png: Path, bins: int = 40) -> None:
    vals = df[col].dropna()
    if vals.empty:
        return
    plt.figure(figsize=(6, 4))
    plt.hist(vals, bins=bins)
    plt.xlabel(col)
    plt.ylabel("Number of samples")
    plt.title(col)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def make_scatter(df: pd.DataFrame, x: str, y: str, out_png: Path, label_outliers: bool = True) -> None:
    tmp = df[[x, y, "sample_id"]].dropna().copy()
    if tmp.empty:
        return

    plt.figure(figsize=(6, 5))
    plt.scatter(tmp[x], tmp[y], s=18)

    if label_outliers and len(tmp) > 0:
        zx = robust_z(tmp[x]).abs()
        zy = robust_z(tmp[y]).abs()
        label_mask = (zx > 3) | (zy > 3)
        for _, row in tmp[label_mask].iterrows():
            plt.annotate(row["sample_id"], (row[x], row[y]), fontsize=7)

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{y} vs {x}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def main():
    if not SUMMARY_FILE.exists():
        raise FileNotFoundError(f"Missing summary file: {SUMMARY_FILE}")

    df = pd.read_csv(SUMMARY_FILE, sep="\t")
    if "sample_id" not in df.columns:
        raise ValueError("Expected a sample_id column in summary table")

    # Keep only columns that exist
    metrics = [c for c in KEY_METRICS if c in df.columns]

    # Robust z-scores for main metrics
    for col in metrics:
        df[f"{col}__robust_z"] = robust_z(df[col])

    # Flagging logic
    # Low-confidence / suspiciously noisy samples
    add_flag(df, "flag_low_pass_count", df["n_pass_total"] < 30 if "n_pass_total" in df else False, "flag_low_pass_count")
    add_flag(df, "flag_very_low_pass_count", df["n_pass_total"] < 15 if "n_pass_total" in df else False, "flag_very_low_pass_count")

    add_flag(
        df,
        "flag_low_nonsyn_rare",
        df["n_pass_nonsyn_rare_1e_4"] < 10 if "n_pass_nonsyn_rare_1e_4" in df else False,
        "flag_low_nonsyn_rare",
    )

    add_flag(
        df,
        "flag_low_frac_pass",
        df["frac_pass"] < 0.10 if "frac_pass" in df else False,
        "flag_low_frac_pass",
    )

    add_flag(
        df,
        "flag_many_common_variants",
        df["frac_pass_with_max_af_gt_0_01"] > 0.10 if "frac_pass_with_max_af_gt_0_01" in df else False,
        "flag_many_common_variants",
    )

    add_flag(
        df,
        "flag_many_low_vaf",
        df["frac_vaf_ge_0_25"] < 0.20 if "frac_vaf_ge_0_25" in df else False,
        "flag_many_low_vaf",
    )

    add_flag(
        df,
        "flag_low_depth",
        df["dp_median_pass"] < 100 if "dp_median_pass" in df else False,
        "flag_low_depth",
    )

    add_flag(
        df,
        "flag_extreme_indel_fraction",
        (
            (df["n_pass_total"] > 0) &
            ((df["n_pass_indel"] / df["n_pass_total"]) > 0.35)
        ) if {"n_pass_indel", "n_pass_total"}.issubset(df.columns) else False,
        "flag_extreme_indel_fraction",
    )

    # High-burden outliers
    if "burden_nonsyn_rare_1e_4_perMb" in df.columns:
        z = df["burden_nonsyn_rare_1e_4_perMb__robust_z"].abs()
        add_flag(df, "flag_high_rare_burden_outlier", z > 3.5, "flag_high_rare_burden_outlier")

    if "n_pass_total" in df.columns:
        z = df["n_pass_total__robust_z"].abs()
        add_flag(df, "flag_pass_count_outlier", z > 3.5, "flag_pass_count_outlier")

    if "n_pass_indel" in df.columns and "n_pass_total" in df.columns:
        df["pass_indel_fraction"] = np.where(df["n_pass_total"] > 0, df["n_pass_indel"] / df["n_pass_total"], np.nan)
        df["pass_indel_fraction__robust_z"] = robust_z(df["pass_indel_fraction"])
    else:
        df["pass_indel_fraction"] = np.nan
        df["pass_indel_fraction__robust_z"] = np.nan

    # Aggregate flag count
    flag_cols = [c for c in df.columns if c.startswith("flag_")]
    df["n_qc_flags"] = df[flag_cols].sum(axis=1)

    # Suggested status
    df["qc_status"] = "PASS"
    df.loc[df["n_qc_flags"] >= 1, "qc_status"] = "REVIEW"
    df.loc[df["n_qc_flags"] >= 3, "qc_status"] = "FAIL_CANDIDATE"

    # Human-readable reason string
    df["qc_reasons"] = df[flag_cols].apply(
        lambda row: ",".join([c for c in flag_cols if row[c] == 1]),
        axis=1
    )

    # Write QC table
    df.to_csv(QC_TABLE, sep="\t", index=False)

    # Build text report
    lines = []
    lines.append("VEP summary QC report")
    lines.append("=" * 80)
    lines.append(f"Input file: {SUMMARY_FILE}")
    lines.append(f"N samples: {df.shape[0]}")
    lines.append("")

    lines.append("Key metric summaries")
    lines.append("-" * 80)
    for col in metrics + ["pass_indel_fraction"]:
        if col not in df.columns:
            continue
        ss = summarize_series(df[col])
        lines.append(
            f"{col:30s} "
            f"n={ss['n']:3d} "
            f"mean={ss['mean']:.3f} "
            f"median={ss['median']:.3f} "
            f"p05={ss['p05']:.3f} "
            f"p95={ss['p95']:.3f} "
            f"min={ss['min']:.3f} "
            f"max={ss['max']:.3f}"
        )

    lines.append("")
    lines.append("QC flag counts")
    lines.append("-" * 80)
    for col in flag_cols:
        lines.append(f"{col:30s} {int(df[col].sum())}")

    lines.append("")
    lines.append("QC status counts")
    lines.append("-" * 80)
    for status, n in df["qc_status"].value_counts(dropna=False).items():
        lines.append(f"{status:20s} {int(n)}")

    lines.append("")
    lines.append("Top flagged samples")
    lines.append("-" * 80)
    flagged = df.sort_values(["n_qc_flags", "n_pass_total"], ascending=[False, True]).copy()
    flagged = flagged.loc[flagged["n_qc_flags"] > 0, ["sample_id", "qc_status", "n_qc_flags", "qc_reasons",
                                                      "n_pass_total", "n_pass_nonsyn_rare_1e_4",
                                                      "frac_pass", "dp_median_pass", "vaf_median_pass",
                                                      "frac_pass_with_max_af_gt_0_01"]]
    if flagged.empty:
        lines.append("No flagged samples.")
    else:
        for _, row in flagged.head(40).iterrows():
            lines.append(
                f"{row['sample_id']:15s} "
                f"{row['qc_status']:14s} "
                f"flags={int(row['n_qc_flags'])} "
                f"n_pass={row.get('n_pass_total', np.nan)} "
                f"rare_nonsyn={row.get('n_pass_nonsyn_rare_1e_4', np.nan)} "
                f"frac_pass={row.get('frac_pass', np.nan):.3f} "
                f"dp_med={row.get('dp_median_pass', np.nan):.1f} "
                f"vaf_med={row.get('vaf_median_pass', np.nan):.3f} "
                f"common_frac={row.get('frac_pass_with_max_af_gt_0_01', np.nan):.3f} "
                f"reasons={row['qc_reasons']}"
            )

    with open(QC_REPORT, "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines[:60]))
    print(f"\nWrote QC table:   {QC_TABLE}")
    print(f"Wrote QC report:  {QC_REPORT}")

    # Plots
    plot_cols = [
        "n_pass_total",
        "n_pass_nonsyn",
        "n_pass_nonsyn_rare_1e_4",
        "burden_nonsyn_rare_1e_4_perMb",
        "frac_pass",
        "vaf_median_pass",
        "dp_median_pass",
        "frac_pass_with_max_af_gt_0_01",
        "pass_indel_fraction",
    ]

    for col in plot_cols:
        if col in df.columns:
            make_hist(df, col, OUT_DIR / f"{col}.hist.png")

    scatter_pairs = [
        ("n_pass_total", "n_pass_nonsyn_rare_1e_4"),
        ("dp_median_pass", "vaf_median_pass"),
        ("frac_pass_with_max_af_gt_0_01", "n_pass_nonsyn_rare_1e_4"),
        ("pass_indel_fraction", "n_pass_total"),
        ("burden_nonsyn_rare_1e_4_perMb", "frac_pass_with_max_af_gt_0_01"),
    ]

    for x, y in scatter_pairs:
        if x in df.columns and y in df.columns:
            make_scatter(df, x, y, OUT_DIR / f"{y}_vs_{x}.png")


if __name__ == "__main__":
    main()
