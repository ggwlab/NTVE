#!/usr/bin/env python3
"""Run the direct harmonized workflow: STAR alignment and featureCounts.

Workflow reconstructed from:
- star_on_harmonized_data.ipynb (STAR on harmonized index)
- run_harmonized_submission_direct.ipynb (featureCounts, s=0 and s=2)
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

PAIR_PATTERNS = [
    re.compile(r"(.+)_R1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
    re.compile(r"(.+)_1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
    re.compile(r"(.+)\.R1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("--mode", choices=["star", "featurecounts", "all"], default="all")

    p.add_argument("--workflow-root", type=Path, default=Path("."), help="Base directory for outputs")
    p.add_argument("--project", help="Project label used in STAR output path (e.g., NTVE)")
    p.add_argument("--species", default="harmonized", help="Species label used in STAR output path")

    p.add_argument("--raw-data-root", type=Path, help="Root directory to scan for FASTQ(.gz) files")
    p.add_argument("--star-bin", default="STAR", help="Path or command name for STAR")
    p.add_argument("--star-index", type=Path, help="STAR index directory")
    p.add_argument("--star-threads", type=int, default=8)
    p.add_argument("--executor", choices=["local", "slurm"], default="local", help="How STAR commands are executed")
    p.add_argument("--slurm-options", default="--time=03:00:00 --mem=110G -n 8 --job-name STAR_RNA")
    p.add_argument("--wait-for-star", action="store_true", help="Wait for submitted STAR SLURM jobs to finish")
    p.add_argument("--poll-seconds", type=int, default=30)

    p.add_argument("--featurecounts-bin", default="featureCounts", help="Path or command name for featureCounts")
    p.add_argument("--annotation-gtf", type=Path, help="GTF annotation file for featureCounts")
    p.add_argument("--bam-glob", help="Optional BAM glob override for featureCounts")
    p.add_argument("--featurecounts-output-dir", type=Path, default=Path("featurecounts_results"))
    p.add_argument("--featurecounts-threads", type=int, default=32)
    p.add_argument("--feature-type", default="exon")
    p.add_argument("--group-by", default="gene_id")

    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def binary_exists(cmd: str) -> bool:
    return Path(cmd).exists() or shutil.which(cmd) is not None


def q(value: str | Path) -> str:
    return shlex.quote(str(value))


def sample_name_from_r1(filename: str) -> str | None:
    for rx in PAIR_PATTERNS:
        m = rx.match(filename)
        if m:
            return m.group(1)
    return None


def derive_r2_filename(r1_name: str) -> str | None:
    if "_R1" in r1_name:
        return r1_name.replace("_R1", "_R2", 1)
    if "_1" in r1_name:
        return r1_name.replace("_1", "_2", 1)
    if ".R1" in r1_name:
        return r1_name.replace(".R1", ".R2", 1)
    return None


def discover_fastq_pairs(raw_root: Path) -> tuple[dict[str, tuple[Path, Path]], list[str]]:
    fastq_files = sorted(raw_root.glob("**/*.fastq.gz")) + sorted(raw_root.glob("**/*.fq.gz"))
    by_name = {p.name: p for p in fastq_files}

    pairs: dict[str, tuple[Path, Path]] = {}
    unmatched_r1: list[str] = []

    for p in fastq_files:
        sample = sample_name_from_r1(p.name)
        if sample is None:
            continue
        r2_name = derive_r2_filename(p.name)
        if r2_name and r2_name in by_name:
            pairs[sample] = (p, by_name[r2_name])
        else:
            unmatched_r1.append(p.name)

    return pairs, unmatched_r1


def star_output_dir(workflow_root: Path, project: str, species: str) -> Path:
    return workflow_root / "star_alignments" / f"Project_{project}" / f"star_{species}"


def build_star_command(
    star_bin: str,
    star_index: Path,
    r1: Path,
    r2: Path,
    output_prefix: Path,
    threads: int,
) -> list[str]:
    return [
        star_bin,
        "--readFilesCommand",
        "zcat",
        "--genomeDir",
        str(star_index),
        "--readFilesIn",
        str(r1),
        str(r2),
        "--runThreadN",
        str(threads),
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        "--outFileNamePrefix",
        str(output_prefix),
        "--quantMode",
        "GeneCounts",
        "TranscriptomeSAM",
    ]


def submit_sbatch(command: list[str], slurm_options: str, sample_name: str) -> str | None:
    cmd_string = " ".join(q(x) for x in command)

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".sh") as tmp:
        tmp.write("#!/bin/bash\n")
        tmp.write(f"#SBATCH {slurm_options}\n")
        tmp.write(f"# sample={sample_name}\n")
        tmp.write("set -euo pipefail\n")
        tmp.write(f"{cmd_string}\n")
        tmp_path = Path(tmp.name)

    try:
        result = subprocess.run(["sbatch", str(tmp_path)], capture_output=True, text=True, check=False)
        if result.returncode != 0:
            return None
        return result.stdout.strip().split()[-1]
    finally:
        tmp_path.unlink(missing_ok=True)


def wait_for_slurm_jobs(job_ids: list[str], poll_seconds: int) -> None:
    pending = set(job_ids)
    while pending:
        query = subprocess.run(["squeue", "-h", "-j", ",".join(sorted(pending)), "-o", "%A"], capture_output=True, text=True, check=False)
        if query.returncode != 0:
            raise RuntimeError(f"squeue failed: {query.stderr.strip()}")

        alive = {line.strip() for line in query.stdout.splitlines() if line.strip()}
        finished = pending - alive
        for job in sorted(finished):
            print(f"STAR finished: job {job}")
        pending = alive

        if pending:
            print(f"Waiting for {len(pending)} STAR job(s)...")
            time.sleep(poll_seconds)


def run_star_stage(args: argparse.Namespace) -> dict:
    if args.raw_data_root is None:
        raise ValueError("--raw-data-root is required for mode 'star' or 'all'")
    if args.star_index is None:
        raise ValueError("--star-index is required for mode 'star' or 'all'")
    if not args.project:
        raise ValueError("--project is required for mode 'star' or 'all'")
    if not args.raw_data_root.exists():
        raise FileNotFoundError(f"Raw data root not found: {args.raw_data_root}")
    if not args.star_index.exists():
        raise FileNotFoundError(f"STAR index not found: {args.star_index}")
    if not binary_exists(args.star_bin):
        raise FileNotFoundError(f"STAR not found: {args.star_bin}")

    if args.executor == "slurm":
        for tool in ["sbatch"] + (["squeue"] if args.wait_for_star else []):
            if not binary_exists(tool):
                raise FileNotFoundError(f"Required SLURM tool not found: {tool}")

    pairs, unmatched_r1 = discover_fastq_pairs(args.raw_data_root)
    print(f"STAR sample pairs discovered: {len(pairs)}")
    print(f"Unmatched R1 files: {len(unmatched_r1)}")

    out_dir = star_output_dir(args.workflow_root, args.project, args.species)
    out_dir.mkdir(parents=True, exist_ok=True)

    submitted_jobs: list[tuple[str, str]] = []
    failed_jobs: list[str] = []
    finished_samples: list[str] = []
    star_commands: dict[str, str] = {}

    for sample_name, (r1, r2) in sorted(pairs.items()):
        output_prefix = out_dir / f"{sample_name}_untrimmed_{args.species}_"
        cmd = build_star_command(args.star_bin, args.star_index, r1, r2, output_prefix, args.star_threads)
        star_commands[sample_name] = " ".join(q(x) for x in cmd)

        if args.dry_run:
            print(f"DRY RUN STAR {sample_name}")
            continue

        if args.executor == "local":
            subprocess.run(cmd, check=True)
            finished_samples.append(sample_name)
            print(f"STAR completed: {sample_name}")
        else:
            job_id = submit_sbatch(cmd, args.slurm_options, sample_name)
            if job_id:
                submitted_jobs.append((sample_name, job_id))
                print(f"STAR submitted: {sample_name} -> {job_id}")
            else:
                failed_jobs.append(sample_name)
                print(f"STAR failed to submit: {sample_name}")

    if args.executor == "slurm" and args.wait_for_star and submitted_jobs and not args.dry_run:
        wait_for_slurm_jobs([job_id for _, job_id in submitted_jobs], args.poll_seconds)

    return {
        "star_pairs_discovered": len(pairs),
        "unmatched_r1_files": unmatched_r1,
        "star_output_dir": str(out_dir),
        "star_submitted_jobs": submitted_jobs,
        "star_failed_jobs": failed_jobs,
        "star_finished_samples": finished_samples,
        "star_commands": star_commands,
    }


def discover_bams_for_featurecounts(args: argparse.Namespace) -> list[str]:
    if args.bam_glob:
        bam_files = sorted(glob.glob(args.bam_glob))
        if not bam_files:
            raise FileNotFoundError(f"No BAM files matched glob: {args.bam_glob}")
        return bam_files

    if not args.project:
        raise ValueError("--project is required for featureCounts when --bam-glob is not provided")

    pattern = str(star_output_dir(args.workflow_root, args.project, args.species) / "*.bam")
    bam_files = sorted(glob.glob(pattern))
    if not bam_files:
        raise FileNotFoundError(f"No BAM files found with default pattern: {pattern}")
    return bam_files


def run_featurecounts_stage(args: argparse.Namespace) -> dict:
    if args.annotation_gtf is None:
        raise ValueError("--annotation-gtf is required for mode 'featurecounts' or 'all'")
    if not args.annotation_gtf.exists():
        raise FileNotFoundError(f"Annotation GTF not found: {args.annotation_gtf}")
    if not binary_exists(args.featurecounts_bin):
        raise FileNotFoundError(f"featureCounts not found: {args.featurecounts_bin}")

    bam_files = discover_bams_for_featurecounts(args)
    print(f"featureCounts BAM files discovered: {len(bam_files)}")

    out_dir = args.workflow_root / args.featurecounts_output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    unstranded_out = out_dir / "harmonized_harmonized_gene_counts_unstranded.txt"
    rv_stranded_out = out_dir / "harmonized_harmonized_gene_counts_rv_stranded.txt"

    common = [
        args.featurecounts_bin,
        "-a",
        str(args.annotation_gtf),
        "-g",
        args.group_by,
        "-t",
        args.feature_type,
        "-p",
        "-T",
        str(args.featurecounts_threads),
    ]

    cmd_unstranded = common + ["-s", "0", "-o", str(unstranded_out)] + bam_files
    cmd_rv = common + ["-s", "2", "-o", str(rv_stranded_out)] + bam_files

    if args.dry_run:
        print("DRY RUN featureCounts")
    else:
        subprocess.run(cmd_unstranded, check=True)
        subprocess.run(cmd_rv, check=True)
        print("featureCounts completed")

    return {
        "featurecounts_bam_count": len(bam_files),
        "featurecounts_bam_files": bam_files,
        "featurecounts_outputs": {
            "unstranded": str(unstranded_out),
            "reverse_stranded": str(rv_stranded_out),
        },
        "featurecounts_commands": {
            "unstranded": " ".join(q(x) for x in cmd_unstranded),
            "reverse_stranded": " ".join(q(x) for x in cmd_rv),
        },
    }


def main() -> int:
    args = parse_args()

    tracking: dict = {
        "mode": args.mode,
        "dry_run": args.dry_run,
        "workflow_root": str(args.workflow_root),
        "project": args.project,
        "species": args.species,
        "executor": args.executor,
        "star_threads": args.star_threads,
        "featurecounts_threads": args.featurecounts_threads,
    }

    if args.mode in {"star", "all"}:
        tracking.update(run_star_stage(args))

    if args.mode in {"featurecounts", "all"}:
        if args.mode == "all" and args.executor == "slurm" and not args.wait_for_star and not args.dry_run:
            raise RuntimeError(
                "mode=all with executor=slurm requires --wait-for-star, "
                "otherwise featureCounts may start before STAR outputs exist"
            )
        tracking.update(run_featurecounts_stage(args))

    tracking_dir = args.workflow_root
    tracking_dir.mkdir(parents=True, exist_ok=True)
    tracking_file = tracking_dir / "direct_workflow_tracking.json"
    tracking_file.write_text(json.dumps(tracking, indent=2))

    print(f"Tracking file: {tracking_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
