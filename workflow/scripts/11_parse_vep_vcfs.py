#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from pathlib import Path


GENES_OF_INTEREST = ["TP53", "FLT3", "NRAS", "KRAS"]
OPTIONAL_DRIVER_GENES = [
    "DNMT3A", "NPM1", "IDH1", "IDH2", "RUNX1", "TET2", "ASXL1",
    "CEBPA", "KIT", "WT1", "PTPN11", "SRSF2", "SF3B1", "U2AF1",
    "STAG2", "BCOR", "EZH2", "PHF6", "CBL", "KRAS", "NRAS", "FLT3", "TP53"
]
ALL_FLAG_GENES = sorted(set(GENES_OF_INTEREST + OPTIONAL_DRIVER_GENES))

NONSYN_TERMS = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
}

TRUNC_TERMS = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

SPLICE_TERMS = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
}

SYN_TERMS = {
    "synonymous_variant",
}

RARE_THRESHOLDS = {
    "1e-3": 1e-3,
    "1e-4": 1e-4,
}


def safe_float(x):
    try:
        if x is None or x == "" or x == ".":
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def parse_panel_mb(panel_size_file: Path) -> float:
    if not panel_size_file.exists():
        raise FileNotFoundError(f"Panel size file not found: {panel_size_file}")

    panel_mb = None
    with open(panel_size_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            key, value = parts
            if key == "panel_mb":
                panel_mb = float(value)
                break

    if panel_mb is None:
        raise ValueError(f"Could not parse panel_mb from {panel_size_file}")

    return panel_mb


def parse_vep_csq_header(vcf: VCF) -> list[str]:
    hdr = vcf.raw_header.splitlines()
    for line in hdr:
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r'Format: (.+?)">', line)
            if m:
                return m.group(1).split("|")
    raise RuntimeError("Could not find CSQ header format in VCF")


def extract_best_csq(csq_str: str | None, csq_fields: list[str]) -> dict[str, str]:
    if csq_str is None:
        return {}
    first = csq_str.split(",")[0]
    parts = first.split("|")
    if len(parts) < len(csq_fields):
        parts += [""] * (len(csq_fields) - len(parts))
    return dict(zip(csq_fields, parts))


def get_filter_label(var) -> str:
    if var.FILTER is None or var.FILTER == "PASS":
        return "PASS"
    return str(var.FILTER)


def extract_dp(var) -> float:
    try:
        dp = var.format("DP")
        if dp is not None and dp.size > 0:
            if dp.ndim == 2:
                return safe_float(dp[0][0])
            return safe_float(dp[0])
    except Exception:
        pass

    try:
        return safe_float(var.INFO.get("DP"))
    except Exception:
        return np.nan


def extract_ad(var) -> tuple[float, float]:
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


def extract_af(var) -> float:
    try:
        af = var.format("AF")
        if af is not None and af.size > 0:
            if af.ndim == 2:
                return safe_float(af[0][0])
            return safe_float(af[0])
    except Exception:
        pass
    return np.nan


def compute_vaf(var) -> tuple[float, float, float, float]:
    dp = extract_dp(var)
    ad_ref, ad_alt = extract_ad(var)
    af = extract_af(var)
    vaf = af
    if np.isnan(vaf) and not np.isnan(ad_ref) and not np.isnan(ad_alt):
        denom = ad_ref + ad_alt
        if denom > 0:
            vaf = ad_alt / denom
    return dp, ad_ref, ad_alt, vaf


def classify_variant(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "INDEL"


def parse_consequence_set(consequence: str) -> set[str]:
    if not consequence:
        return set()
    terms = set()
    for block in consequence.split("&"):
        block = block.strip()
        if block:
            terms.add(block)
    return terms


def consequence_flags(consequence: str) -> dict[str, bool]:
    terms = parse_consequence_set(consequence)
    return {
        "is_syn": len(terms & SYN_TERMS) > 0,
        "is_nonsyn": len(terms & NONSYN_TERMS) > 0,
        "is_trunc": len(terms & TRUNC_TERMS) > 0,
        "is_splice": len(terms & SPLICE_TERMS) > 0,
    }


def parse_existing_variation(existing_variation: str) -> dict[str, object]:
    if not existing_variation:
        return {
            "existing_variation": "",
            "has_existing_variation": False,
            "n_existing_ids": 0,
            "has_rsid": False,
            "has_cosmic": False,
        }

    ids = [x for x in existing_variation.split("&") if x]
    return {
        "existing_variation": existing_variation,
        "has_existing_variation": len(ids) > 0,
        "n_existing_ids": len(ids),
        "has_rsid": any(x.startswith("rs") for x in ids),
        "has_cosmic": any(x.startswith("COS") for x in ids),
    }


def is_rare(max_af: float, threshold: float) -> bool:
    return np.isnan(max_af) or (max_af < threshold)


def nanmedian_safe(values):
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.nanmedian(arr))


def nanmax_safe(values):
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.nanmax(arr))


def nanmean_safe(values):
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.nanmean(arr))


def naniqr_safe(values):
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.nanpercentile(arr, 75) - np.nanpercentile(arr, 25))


def summarize_sample(sample_id, panel_mb, rows, gene_hits, filter_counter):
    if len(rows) == 0:
        return {"sample_id": sample_id, "panel_mb": panel_mb}

    df = pd.DataFrame(rows)

    n_total = len(df)
    n_total_snv = int((df["variant_type"] == "SNV").sum())
    n_total_indel = int((df["variant_type"] == "INDEL").sum())

    df_pass = df[df["filter"] == "PASS"].copy()
    n_pass_total = len(df_pass)
    n_pass_snv = int((df_pass["variant_type"] == "SNV").sum())
    n_pass_indel = int((df_pass["variant_type"] == "INDEL").sum())

    n_pass_syn = int(df_pass["is_syn"].sum())
    n_pass_nonsyn = int(df_pass["is_nonsyn"].sum())
    n_pass_trunc = int(df_pass["is_trunc"].sum())
    n_pass_splice = int(df_pass["is_splice"].sum())

    rare_counts = {}
    for label, thr in RARE_THRESHOLDS.items():
        rare_counts[f"n_pass_nonsyn_rare_{label.replace('-', '_')}"] = int(
            ((df_pass["is_nonsyn"]) & ((df_pass["max_af"].isna()) | (df_pass["max_af"] < thr))).sum()
        )

    pass_vaf = df_pass["vaf"].dropna().tolist()
    pass_dp = df_pass["dp"].dropna().tolist()
    pass_ad_alt = df_pass["ad_alt"].dropna().tolist()

    summary = {
        "sample_id": sample_id,
        "panel_mb": panel_mb,
        "n_total": n_total,
        "n_total_snv": n_total_snv,
        "n_total_indel": n_total_indel,
        "n_pass_total": n_pass_total,
        "n_pass_snv": n_pass_snv,
        "n_pass_indel": n_pass_indel,
        "n_pass_syn": n_pass_syn,
        "n_pass_nonsyn": n_pass_nonsyn,
        "n_pass_trunc": n_pass_trunc,
        "n_pass_splice": n_pass_splice,
        **rare_counts,
        "burden_total_perMb": n_pass_total / panel_mb,
        "burden_nonsyn_perMb": n_pass_nonsyn / panel_mb,
        "burden_nonsyn_rare_1e_3_perMb": rare_counts["n_pass_nonsyn_rare_1e_3"] / panel_mb,
        "burden_nonsyn_rare_1e_4_perMb": rare_counts["n_pass_nonsyn_rare_1e_4"] / panel_mb,
        "frac_pass": (n_pass_total / n_total) if n_total > 0 else np.nan,
        "vaf_median_pass": nanmedian_safe(pass_vaf),
        "vaf_max_pass": nanmax_safe(pass_vaf),
        "vaf_iqr_pass": naniqr_safe(pass_vaf),
        "frac_vaf_ge_0_25": nanmean_safe((df_pass["vaf"] >= 0.25).astype(float).tolist()) if n_pass_total > 0 else np.nan,
        "frac_vaf_ge_0_40": nanmean_safe((df_pass["vaf"] >= 0.40).astype(float).tolist()) if n_pass_total > 0 else np.nan,
        "n_low_vaf_lt_0_05": int((df_pass["vaf"] < 0.05).sum()),
        "dp_median_pass": nanmedian_safe(pass_dp),
        "dp_iqr_pass": naniqr_safe(pass_dp),
        "ad_alt_median_pass": nanmedian_safe(pass_ad_alt),
        "frac_pass_with_existing_variation": nanmean_safe(df_pass["has_existing_variation"].astype(float).tolist()) if n_pass_total > 0 else np.nan,
        "frac_pass_with_rsid": nanmean_safe(df_pass["has_rsid"].astype(float).tolist()) if n_pass_total > 0 else np.nan,
        "frac_pass_with_max_af": nanmean_safe(df_pass["max_af"].notna().astype(float).tolist()) if n_pass_total > 0 else np.nan,
        "frac_pass_with_max_af_gt_0_01": nanmean_safe((df_pass["max_af"] > 0.01).astype(float).tolist()) if n_pass_total > 0 else np.nan,
    }

    summary["n_filter_PASS"] = filter_counter.get("PASS", 0)
    for key, val in sorted(filter_counter.items()):
        safe_key = re.sub(r"[^A-Za-z0-9_]+", "_", key)
        summary[f"n_filter_{safe_key}"] = val

    for gene in ALL_FLAG_GENES:
        summary[f"{gene}_mut"] = int(gene_hits.get(gene, 0))

    summary["RAS_mut"] = int(summary.get("NRAS_mut", 0) or summary.get("KRAS_mut", 0))
    return summary


def main():
    parser = argparse.ArgumentParser(description="Parse per-sample VEP VCFs into long and summary tables.")
    parser.add_argument("--vep-dir", required=True)
    parser.add_argument("--panel-size-file", required=True)
    parser.add_argument("--output-long", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    vep_dir = Path(args.vep_dir)
    panel_size_file = Path(args.panel_size_file)
    output_long = Path(args.output_long)
    output_summary = Path(args.output_summary)

    if not vep_dir.exists():
        raise FileNotFoundError(f"VEP directory not found: {vep_dir}")

    output_long.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)

    panel_mb = parse_panel_mb(panel_size_file)
    vcf_files = sorted(vep_dir.glob("*.filtered.vep.vcf.gz"))

    if not vcf_files:
        raise RuntimeError(f"No VEP VCFs found in {vep_dir}")

    all_rows = []
    all_summaries = []

    for i, vcf_path in enumerate(vcf_files, start=1):
        sample_id = vcf_path.name.replace(".filtered.vep.vcf.gz", "")
        print(f"[{i}/{len(vcf_files)}] Parsing {sample_id}")

        vcf = VCF(str(vcf_path))
        csq_fields = parse_vep_csq_header(vcf)

        sample_rows = []
        gene_hits = {g: 0 for g in ALL_FLAG_GENES}
        filter_counter = Counter()

        for var in vcf:
            filter_label = get_filter_label(var)
            filter_counter[filter_label] += 1

            alt_list = list(var.ALT)
            alt = alt_list[0] if alt_list else "."

            dp, ad_ref, ad_alt, vaf = compute_vaf(var)
            variant_type = classify_variant(var.REF, alt)

            csq = extract_best_csq(var.INFO.get("CSQ"), csq_fields)

            consequence = csq.get("Consequence", "")
            symbol = csq.get("SYMBOL", "")
            existing_variation = csq.get("Existing_variation", "")
            max_af = safe_float(csq.get("MAX_AF", ""))

            flags = consequence_flags(consequence)
            existing_info = parse_existing_variation(existing_variation)

            if filter_label == "PASS" and flags["is_nonsyn"] and symbol in gene_hits:
                gene_hits[symbol] = 1

            row = {
                "sample_id": sample_id,
                "chrom": var.CHROM,
                "pos": var.POS,
                "ref": var.REF,
                "alt": alt,
                "variant_type": variant_type,
                "filter": filter_label,
                "dp": safe_float(dp),
                "ad_ref": safe_float(ad_ref),
                "ad_alt": safe_float(ad_alt),
                "vaf": safe_float(vaf),
                "gene": symbol,
                "gene_id": csq.get("Gene", ""),
                "feature_type": csq.get("Feature_type", ""),
                "feature": csq.get("Feature", ""),
                "biotype": csq.get("BIOTYPE", ""),
                "consequence": consequence,
                "impact": csq.get("IMPACT", ""),
                "variant_class": csq.get("VARIANT_CLASS", ""),
                "canonical": csq.get("CANONICAL", ""),
                "hgvsc": csq.get("HGVSc", ""),
                "hgvsp": csq.get("HGVSp", ""),
                "sift": csq.get("SIFT", ""),
                "polyphen": csq.get("PolyPhen", ""),
                "clin_sig": csq.get("CLIN_SIG", ""),
                "somatic": csq.get("SOMATIC", ""),
                **existing_info,
                "af_1kg": safe_float(csq.get("AF", "")),
                "gnomade_af": safe_float(csq.get("gnomADe_AF", "")),
                "gnomadg_af": safe_float(csq.get("gnomADg_AF", "")),
                "max_af": max_af,
                "max_af_pops": csq.get("MAX_AF_POPS", ""),
                "is_syn": flags["is_syn"],
                "is_nonsyn": flags["is_nonsyn"],
                "is_trunc": flags["is_trunc"],
                "is_splice": flags["is_splice"],
            }

            for label, thr in RARE_THRESHOLDS.items():
                row[f"is_rare_{label.replace('-', '_')}"] = is_rare(max_af, thr)

            sample_rows.append(row)
            all_rows.append(row)

        summary = summarize_sample(
            sample_id=sample_id,
            panel_mb=panel_mb,
            rows=sample_rows,
            gene_hits=gene_hits,
            filter_counter=filter_counter,
        )
        all_summaries.append(summary)

    df_long = pd.DataFrame(all_rows)
    df_summary = pd.DataFrame(all_summaries).sort_values("sample_id")

    with gzip.open(output_long, "wt") as f:
        df_long.to_csv(f, sep="\t", index=False)

    df_summary.to_csv(output_summary, sep="\t", index=False)

    print(f"Wrote long table: {output_long}")
    print(f"Wrote summary table: {output_summary}")


if __name__ == "__main__":
    main()
