#!/usr/bin/env python3

from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from cyvcf2 import VCF


PROJECT_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq")
INPUT_DIR = PROJECT_DIR / "results" / "mutect2_norm"
OUTPUT_DIR = PROJECT_DIR / "results" / "analysis"
SAMPLE_LIST = PROJECT_DIR / "metadata" / "all_sample_ids.txt"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_LONG = OUTPUT_DIR / "variants_long_pass.tsv.gz"
OUT_SUMMARY = OUTPUT_DIR / "sample_mutect2_summary_raw.tsv"


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan


def is_pass(variant):
    return (variant.FILTER is None) or (variant.FILTER == "PASS")


def extract_dp(var):
    try:
        dp = var.format("DP")
        if dp is not None and dp.size > 0:
            if dp.ndim == 2:
                return safe_float(dp[0][0])
            return safe_float(dp[0])
    except Exception:
        pass
    return np.nan


def extract_ad(var):
    ad_ref, ad_alt = np.nan, np.nan
    try:
        ad = var.format("AD")
        if ad is not None and ad.size > 0:
            vals = ad[0]
            if len(vals) >= 2:
                ad_ref = safe_float(vals[0])
                ad_alt = safe_float(vals[1])
    except Exception:
        pass
    return ad_ref, ad_alt


def extract_af(var):
    try:
        af = var.format("AF")
        if af is not None and af.size > 0:
            if af.ndim == 2:
                return safe_float(af[0][0])
            return safe_float(af[0])
    except Exception:
        pass
    return np.nan


def classify_variant(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "INDEL"


def main():
    rows = []
    summaries = []

    with open(SAMPLE_LIST) as f:
        sample_ids = [x.strip() for x in f if x.strip()]

    for sample_id in sample_ids:
        vcf_path = INPUT_DIR / f"{sample_id}.filtered.norm.vcf.gz"

        if not vcf_path.exists():
            print(f"WARNING: missing {vcf_path}")
            continue

        print(f"Parsing {sample_id} ...")
        vcf = VCF(str(vcf_path))

        n_pass_total = 0
        n_pass_snv = 0
        n_pass_indel = 0

        vafs = []
        dps = []

        for var in vcf:
            if not is_pass(var):
                continue

            alt_list = list(var.ALT)
            if len(alt_list) != 1:
                # after bcftools norm -m -any, this should usually be 1
                # but we guard anyway
                for alt in alt_list:
                    pass

            alt = alt_list[0] if alt_list else "."

            dp = extract_dp(var)
            ad_ref, ad_alt = extract_ad(var)
            af = extract_af(var)

            vaf = af
            if np.isnan(vaf) and not np.isnan(ad_ref) and not np.isnan(ad_alt):
                denom = ad_ref + ad_alt
                if denom > 0:
                    vaf = ad_alt / denom

            var_type = classify_variant(var.REF, alt)

            n_pass_total += 1
            if var_type == "SNV":
                n_pass_snv += 1
            else:
                n_pass_indel += 1

            if not np.isnan(vaf):
                vafs.append(vaf)
            if not np.isnan(dp):
                dps.append(dp)

            rows.append(
                {
                    "sample_id": sample_id,
                    "chrom": var.CHROM,
                    "pos": var.POS,
                    "ref": var.REF,
                    "alt": alt,
                    "variant_type": var_type,
                    "filter": "PASS",
                    "dp": safe_float(dp),
                    "ad_ref": safe_float(ad_ref),
                    "ad_alt": safe_float(ad_alt),
                    "af": safe_float(af),
                    "vaf": safe_float(vaf),
                }
            )

        vaf_arr = np.array(vafs) if vafs else np.array([np.nan])
        dp_arr = np.array(dps) if dps else np.array([np.nan])

        summaries.append(
            {
                "sample_id": sample_id,
                "n_pass_total": n_pass_total,
                "n_pass_snv": n_pass_snv,
                "n_pass_indel": n_pass_indel,
                "vaf_median": float(np.nanmedian(vaf_arr)),
                "vaf_max": float(np.nanmax(vaf_arr)),
                "frac_vaf_ge_0_25": float(np.nanmean(vaf_arr >= 0.25)),
                "dp_median": float(np.nanmedian(dp_arr)),
            }
        )

    df_long = pd.DataFrame(rows)
    df_summary = pd.DataFrame(summaries).sort_values("sample_id")

    with gzip.open(OUT_LONG, "wt") as f:
        df_long.to_csv(f, sep="\t", index=False)

    df_summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Wrote: {OUT_LONG}")
    print(f"Wrote: {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
