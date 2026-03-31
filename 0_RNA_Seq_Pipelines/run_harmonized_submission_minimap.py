#!/usr/bin/env python3
"""Submit harmonized minimap2->filter->salmon jobs to SLURM.

This script is a generalized, path-agnostic version of the original notebook
workflow documented in README.md.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

PAIR_PATTERNS = [
    re.compile(r"(.+)_R1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
    re.compile(r"(.+)_1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
    re.compile(r"(.+)\.R1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
]

INDEX_FILES = {
    "homo_sapiens": "homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz",
    "mus_musculus": "mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz",
    "HS_MM_combined": "homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz",
}


def parse_args() -> argparse.Namespace:
    user = os.environ.get("USER", "user")
    default_scratch = f"/scratch/{user}/rna_harmonized"

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-data-root", type=Path, required=True, help="Root directory to scan for FASTQ(.gz) files")
    p.add_argument("--output-root", type=Path, required=True, help="Output root (salmon outputs + job tracking)")
    p.add_argument("--index-dir", type=Path, required=True, help="Directory containing required transcriptome FASTA(.gz) files")
    p.add_argument("--scratch-dir", type=Path, default=Path(default_scratch), help="Scratch directory for temporary BAMs")

    p.add_argument("--minimap-bin", default="minimap2")
    p.add_argument("--samtools-bin", default="samtools")
    p.add_argument("--salmon-bin", default="salmon")
    p.add_argument("--filter-bin", default="filter_bam_pairs")

    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--complexity-cutoff", type=float, default=0.8)
    p.add_argument("--min-mapped-bases", type=int, default=80)
    p.add_argument("--num-bootstraps", type=int, default=100)
    p.add_argument("--slurm-options", default="--time=20:00:00 --mem=16G -n 8 --job-name harmonized_rna")
    p.add_argument("--dry-run", action="store_true", help="Build jobs but do not submit")
    return p.parse_args()


def binary_exists(cmd: str) -> bool:
    return Path(cmd).exists() or shutil.which(cmd) is not None


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


def discover_fastq_pairs(base_raw_data: Path) -> tuple[dict[str, tuple[Path, Path]], list[str]]:
    fastq_files = sorted(base_raw_data.glob("**/*.fastq.gz")) + sorted(base_raw_data.glob("**/*.fq.gz"))
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


def get_transcriptome_indices(index_dir: Path) -> dict[str, Path]:
    return {name: index_dir / fname for name, fname in INDEX_FILES.items()}


def q(value: str | Path) -> str:
    return shlex.quote(str(value))


def construct_pipeline_command(
    sample_name: str,
    r1_fastq: Path,
    r2_fastq: Path,
    output_root: Path,
    minimap_bin: str,
    samtools_bin: str,
    salmon_bin: str,
    filter_bin: str,
    transcriptome_indices: dict[str, Path],
    complexity_cutoff: float,
    min_mapped_bases: int,
    threads: int,
    num_bootstraps: int,
    scratch_dir: Path,
) -> str:
    steps: list[str] = []

    for index_name, transcriptome_path in transcriptome_indices.items():
        scratch_sample_dir = scratch_dir / index_name / sample_name
        salmon_outdir = output_root / "salmon" / index_name / sample_name / "quant"

        raw_bam = scratch_sample_dir / f"{sample_name}_raw.bam"
        namesorted_bam = scratch_sample_dir / f"{sample_name}_namesorted.bam"
        filtered_bam = scratch_sample_dir / f"{sample_name}_filtered.bam"

        steps.extend(
            [
                f"echo {q(f'Processing {index_name}...')}",
                f"mkdir -p {q(scratch_sample_dir)} {q(salmon_outdir.parent)}",
                (
                    f"{q(minimap_bin)} -ax sr -t {threads} {q(transcriptome_path)} {q(r1_fastq)} {q(r2_fastq)} | "
                    f"{q(samtools_bin)} view -bS -f 0x2 - > {q(raw_bam)}"
                ),
                f"{q(samtools_bin)} sort -n -@ {threads} -o {q(namesorted_bam)} {q(raw_bam)}",
                f"rm -f {q(raw_bam)}",
                (
                    f"{q(filter_bin)} -i {q(namesorted_bam)} -o {q(filtered_bam)} "
                    f"-c {complexity_cutoff} -m {min_mapped_bases}"
                ),
                f"rm -f {q(namesorted_bam)}",
                (
                    f"{q(salmon_bin)} quant -t {q(transcriptome_path)} -l A -a {q(filtered_bam)} "
                    f"-o {q(salmon_outdir)} -p {threads} --numBootstraps {num_bootstraps} --seqBias --gcBias"
                ),
                f"gzip -f {q(salmon_outdir / 'quant.sf')}",
                f"rm -f {q(filtered_bam)}",
                f"rm -rf {q(scratch_sample_dir)}",
            ]
        )

    steps.append(f"echo {q('All indices processed successfully!')}")
    return " && ".join(steps)


def submit_sbatch_job(pipeline_cmd: str, sample_name: str, slurm_options: str) -> str | None:
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".sh") as tmp:
        tmp.write("#!/bin/bash\n")
        tmp.write(f"#SBATCH {slurm_options}\n")
        tmp.write(f"# Sample: {sample_name}\n")
        tmp.write("set -euo pipefail\n")
        tmp.write(f"{pipeline_cmd}\n")
        tmp_path = Path(tmp.name)

    try:
        result = subprocess.run(["sbatch", str(tmp_path)], capture_output=True, text=True, check=False)
        if result.returncode != 0:
            return None
        return result.stdout.strip().split()[-1]
    finally:
        tmp_path.unlink(missing_ok=True)


def main() -> int:
    args = parse_args()

    for required in [args.raw_data_root, args.index_dir]:
        if not required.exists():
            raise FileNotFoundError(f"Required path not found: {required}")

    required_bins = [args.minimap_bin, args.samtools_bin, args.salmon_bin, args.filter_bin]
    missing_bins = [b for b in required_bins if not binary_exists(b)]
    if missing_bins:
        raise FileNotFoundError("Missing required tools: " + ", ".join(missing_bins))

    transcriptome_indices = get_transcriptome_indices(args.index_dir)
    missing_indices = [f"{name} -> {path}" for name, path in transcriptome_indices.items() if not path.exists()]
    if missing_indices:
        raise FileNotFoundError("Missing transcriptome references:\n" + "\n".join(missing_indices))

    sample_pairs, unmatched_r1 = discover_fastq_pairs(args.raw_data_root)
    print(f"Matched FASTQ pairs: {len(sample_pairs)}")
    print(f"Unmatched R1 files: {len(unmatched_r1)}")

    if not sample_pairs:
        print("No FASTQ pairs found. Exiting.")
        return 0

    submitted_jobs: list[tuple[str, str]] = []
    failed_jobs: list[str] = []

    for sample_name, (r1, r2) in sorted(sample_pairs.items()):
        pipeline_cmd = construct_pipeline_command(
            sample_name=sample_name,
            r1_fastq=r1,
            r2_fastq=r2,
            output_root=args.output_root,
            minimap_bin=args.minimap_bin,
            samtools_bin=args.samtools_bin,
            salmon_bin=args.salmon_bin,
            filter_bin=args.filter_bin,
            transcriptome_indices=transcriptome_indices,
            complexity_cutoff=args.complexity_cutoff,
            min_mapped_bases=args.min_mapped_bases,
            threads=args.threads,
            num_bootstraps=args.num_bootstraps,
            scratch_dir=args.scratch_dir,
        )

        if args.dry_run:
            print(f"DRY RUN {sample_name}")
            continue

        job_id = submit_sbatch_job(pipeline_cmd, sample_name, args.slurm_options)
        if job_id:
            submitted_jobs.append((sample_name, job_id))
            print(f"SUBMITTED {sample_name} -> {job_id}")
        else:
            failed_jobs.append(sample_name)
            print(f"FAILED {sample_name}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    tracking = {
        "submitted_jobs": submitted_jobs,
        "failed_jobs": failed_jobs,
        "unmatched_r1_files": unmatched_r1,
        "total_discovered_pairs": len(sample_pairs),
        "total_submitted": len(submitted_jobs),
        "total_failed": len(failed_jobs),
        "raw_data_root": str(args.raw_data_root),
        "output_root": str(args.output_root),
        "index_dir": str(args.index_dir),
        "scratch_dir": str(args.scratch_dir),
        "complexity_cutoff": args.complexity_cutoff,
        "min_mapped_bases": args.min_mapped_bases,
    }

    tracking_file = args.output_root / "job_tracking.json"
    tracking_file.write_text(json.dumps(tracking, indent=2))
    print(f"Tracking file: {tracking_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
