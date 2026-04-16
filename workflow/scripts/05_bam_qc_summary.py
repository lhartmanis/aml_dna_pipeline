#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shlex
import statistics
import subprocess
import time
from datetime import datetime
from pathlib import Path


DEFAULT_SAMPLE_BAM_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/sample_bam")
DEFAULT_QC_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/qc")
DEFAULT_LOG_DIR = Path("/home/leonard.hartmanis/proj/DNA_seq/results/logs")
DEFAULT_INTERVAL_BED = Path("/home/leonard.hartmanis/proj/DNA_seq/intervals/panel_targets_all_hg38.dict_filtered.bed")
DEFAULT_INTERVAL_LIST = Path("/home/leonard.hartmanis/proj/DNA_seq/intervals/panel_targets_all_hg38.interval_list")
DEFAULT_REF_FA = Path("/home/leonard.hartmanis/proj/DNA_seq/ref/GRCh38.p14.genome.fa")


PICARD_FIELDS = [
    "PF_UQ_BASES_ALIGNED",
    "ON_BAIT_BASES",
    "NEAR_BAIT_BASES",
    "OFF_BAIT_BASES",
    "ON_TARGET_BASES",
    "PCT_SELECTED_BASES",
    "PCT_OFF_BAIT",
    "MEAN_BAIT_COVERAGE",
    "MEAN_TARGET_COVERAGE",
    "MEDIAN_TARGET_COVERAGE",
    "FOLD_ENRICHMENT",
    "ZERO_CVG_TARGETS_PCT",
    "PCT_TARGET_BASES_1X",
    "PCT_TARGET_BASES_2X",
    "PCT_TARGET_BASES_10X",
    "PCT_TARGET_BASES_20X",
    "PCT_TARGET_BASES_30X",
    "PCT_TARGET_BASES_50X",
    "PCT_TARGET_BASES_100X",
]


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, logfh) -> None:
    line = f"[{timestamp()}] {msg}"
    print(line)
    logfh.write(line + "\n")
    logfh.flush()


def run_cmd(cmd: list[str], logfh, capture_stdout: bool = False) -> str | None:
    logfh.write("Command:\n")
    logfh.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
    logfh.flush()

    if capture_stdout:
        res = subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=logfh,
        )
        return res.stdout
    else:
        subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh)
        return None

def parse_flagstat(flagstat_text: str) -> dict[str, int]:
    out = {
        "n_total_reads": None,
        "n_mapped_reads": None,
        "n_properly_paired": None,
        "n_duplicate_reads": None,
        "n_paired_in_sequencing": None,
    }

    for line in flagstat_text.strip().splitlines():
        line = line.strip()

        if " in total " in line:
            out["n_total_reads"] = int(line.split()[0])

        elif line.endswith(" mapped (100.00% : N/A)") or " mapped (" in line:
            # THIS is the correct mapped line
            out["n_mapped_reads"] = int(line.split()[0])

        elif " properly paired " in line:
            out["n_properly_paired"] = int(line.split()[0])

        elif " duplicates" in line:
            out["n_duplicate_reads"] = int(line.split()[0])

        elif " paired in sequencing" in line:
            out["n_paired_in_sequencing"] = int(line.split()[0])

    return out


def parse_picard_hsmetrics(metrics_file: Path) -> dict[str, str]:
    with open(metrics_file) as fh:
        lines = [x.rstrip("\n") for x in fh]

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("BAIT_SET\t"):
            header_idx = i
            break

    if header_idx is None or header_idx + 1 >= len(lines):
        raise RuntimeError(f"Could not parse Picard HsMetrics header/row in {metrics_file}")

    headers = lines[header_idx].split("\t")
    values = lines[header_idx + 1].split("\t")
    row = dict(zip(headers, values))

    return {k: row.get(k, "") for k in PICARD_FIELDS}


def summarize_depth(depth_file: Path) -> dict[str, float | None]:
    depths = []
    with open(depth_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                d = int(parts[2])
            except ValueError:
                continue
            depths.append(d)

    if not depths:
        return {
            "mean_target_depth": None,
            "median_target_depth": None,
            "frac_target_bases_1x": None,
            "frac_target_bases_10x": None,
            "frac_target_bases_20x": None,
            "frac_target_bases_50x": None,
            "frac_target_bases_100x": None,
        }

    n = len(depths)
    return {
        "mean_target_depth": sum(depths) / n,
        "median_target_depth": statistics.median(depths),
        "frac_target_bases_1x": sum(d >= 1 for d in depths) / n,
        "frac_target_bases_10x": sum(d >= 10 for d in depths) / n,
        "frac_target_bases_20x": sum(d >= 20 for d in depths) / n,
        "frac_target_bases_50x": sum(d >= 50 for d in depths) / n,
        "frac_target_bases_100x": sum(d >= 100 for d in depths) / n,
    }


def safe_div(num, den):
    if num is None or den in (None, 0):
        return None
    return num / den


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute BAM QC summary and Picard HsMetrics per sample.")
    parser.add_argument("--sample-bam-dir", default=str(DEFAULT_SAMPLE_BAM_DIR))
    parser.add_argument("--qc-dir", default=str(DEFAULT_QC_DIR))
    parser.add_argument("--log-dir", default=str(DEFAULT_LOG_DIR))
    parser.add_argument("--interval-bed", default=str(DEFAULT_INTERVAL_BED))
    parser.add_argument("--interval-list", default=str(DEFAULT_INTERVAL_LIST))
    parser.add_argument("--reference", default=str(DEFAULT_REF_FA))
    parser.add_argument("--sample-id", required=True, help="Run one sample only")
    parser.add_argument("--force", action="store_true", help="Recompute sample QC even if raw outputs already exist")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    sample_bam_dir = Path(args.sample_bam_dir)
    qc_dir = Path(args.qc_dir)
    log_dir = Path(args.log_dir)
    interval_bed = Path(args.interval_bed)
    interval_list = Path(args.interval_list)
    reference = Path(args.reference)
    sample_id = args.sample_id

    if not sample_bam_dir.exists():
        raise FileNotFoundError(f"Sample BAM dir not found: {sample_bam_dir}")
    if not interval_bed.exists():
        raise FileNotFoundError(f"Interval BED not found: {interval_bed}")
    if not interval_list.exists():
        raise FileNotFoundError(f"Interval list not found: {interval_list}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")

    qc_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    bam = sample_bam_dir / f"{sample_id}.markdup.bam"
    bai = sample_bam_dir / f"{sample_id}.markdup.bam.bai"
    sample_log = log_dir / f"{sample_id}.qc.log"

    flagstat_file = qc_dir / f"{sample_id}.flagstat.txt"
    idxstats_file = qc_dir / f"{sample_id}.idxstats.txt"
    depth_file = qc_dir / f"{sample_id}.target_depth.tsv"
    hsmetrics_file = qc_dir / f"{sample_id}.hs_metrics.txt"
    json_file = qc_dir / f"{sample_id}.qc_summary.json"

    if not bam.exists():
        raise FileNotFoundError(f"missing_bam:{bam}")
    if not bai.exists():
        raise FileNotFoundError(f"missing_bai:{bai}")

    if (
        json_file.exists()
        and flagstat_file.exists()
        and idxstats_file.exists()
        and depth_file.exists()
        and hsmetrics_file.exists()
        and not args.force
    ):
        print(f"[{timestamp()}] SKIP {sample_id}: QC outputs already exist")
        return

    if args.dry_run:
        print(f"[{timestamp()}] DRY-RUN {sample_id}")
        return

    start = time.time()

    with open(sample_log, "w") as sample_fh:
        sample_fh.write(f"Start time: {datetime.now()}\n")
        sample_fh.write(f"Sample: {sample_id}\n")
        sample_fh.write(f"BAM: {bam}\n\n")
        sample_fh.flush()

        flagstat_cmd = ["samtools", "flagstat", str(bam)]
        flagstat_text = run_cmd(flagstat_cmd, sample_fh, capture_stdout=True)
        assert flagstat_text is not None
        with open(flagstat_file, "w") as fh:
            fh.write(flagstat_text)
        flagstat = parse_flagstat(flagstat_text)

        idxstats_cmd = ["samtools", "idxstats", str(bam)]
        idxstats_text = run_cmd(idxstats_cmd, sample_fh, capture_stdout=True)
        assert idxstats_text is not None
        with open(idxstats_file, "w") as fh:
            fh.write(idxstats_text)

        # 3332 = 4 (unmapped) + 256 (secondary) + 2048 (supplementary) + 1024 (duplicates)

        on_target_cmd = [
            "samtools", "view",
            "-c",
            "-F", "3332",
            "-L", str(interval_bed),
            str(bam)
        ]

        on_target_text = run_cmd(on_target_cmd, sample_fh, capture_stdout=True)
        assert on_target_text is not None
        n_on_target_reads_rough = int(on_target_text.strip())

        depth_cmd = [
            "samtools", "depth",
            "-a",
            "-b", str(interval_bed),
            str(bam)
        ]
        depth_text = run_cmd(depth_cmd, sample_fh, capture_stdout=True)
        assert depth_text is not None
        with open(depth_file, "w") as fh:
            fh.write(depth_text)
        depth_stats = summarize_depth(depth_file)

        hs_cmd = [
            "picard",
            "-Xmx4g",
            "CollectHsMetrics",
            f"I={bam}",
            f"O={hsmetrics_file}",
            f"R={reference}",
            f"BAIT_INTERVALS={interval_list}",
            f"TARGET_INTERVALS={interval_list}",
            "VALIDATION_STRINGENCY=LENIENT",
            "PER_TARGET_COVERAGE=/dev/null",
        ]
        run_cmd(hs_cmd, sample_fh, capture_stdout=False)
        picard_stats = parse_picard_hsmetrics(hsmetrics_file)

        bam_size_bytes = bam.stat().st_size
        n_total_reads = flagstat["n_total_reads"]
        n_mapped_reads = flagstat["n_mapped_reads"]
        n_properly_paired = flagstat["n_properly_paired"]
        n_duplicate_reads = flagstat["n_duplicate_reads"]
        n_paired_in_seq = flagstat["n_paired_in_sequencing"]

        row = {
            "sample_id": sample_id,
            "bam_path": str(bam),
            "bam_size_bytes": bam_size_bytes,
            "n_total_reads": n_total_reads,
            "n_mapped_reads": n_mapped_reads,
            "frac_mapped": safe_div(n_mapped_reads, n_total_reads),
            "n_properly_paired": n_properly_paired,
            "frac_properly_paired": safe_div(n_properly_paired, n_paired_in_seq),
            "n_duplicate_reads": n_duplicate_reads,
            "frac_duplicates": safe_div(n_duplicate_reads, n_total_reads),
            "n_on_target_reads_rough": n_on_target_reads_rough,
            "frac_on_target_rough": safe_div(n_on_target_reads_rough, n_mapped_reads),
            "mean_target_depth": depth_stats["mean_target_depth"],
            "median_target_depth": depth_stats["median_target_depth"],
            "frac_target_bases_1x": depth_stats["frac_target_bases_1x"],
            "frac_target_bases_10x": depth_stats["frac_target_bases_10x"],
            "frac_target_bases_20x": depth_stats["frac_target_bases_20x"],
            "frac_target_bases_50x": depth_stats["frac_target_bases_50x"],
            "frac_target_bases_100x": depth_stats["frac_target_bases_100x"],
        }
        row.update(picard_stats)

        with open(json_file, "w") as fh:
            json.dump(row, fh, indent=2)

        elapsed = time.time() - start
        sample_fh.write(f"\nJSON summary: {json_file}\n")
        sample_fh.write(f"End time: {datetime.now()}\n")
        sample_fh.write(f"Elapsed seconds: {elapsed:.2f}\n")
        sample_fh.flush()

    print(f"[{timestamp()}] DONE {sample_id}: elapsed={elapsed:.2f}s")


if __name__ == "__main__":
    main()
