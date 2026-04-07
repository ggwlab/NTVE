"""
Supplementary Figure 20a — FastQC/BLAST quality comparison.

Refactors the reproducible part of:
- Suppl20/fastqc_based_analysis.ipynb

This script intentionally uses the cached FastQC-derived FASTA files and cached
BLAST text outputs already present in the repo. It does not submit live BLAST
queries. The goal is to reproduce the notebook's comparison figures and the
intermediate summary table offline.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent


def first_existing_dir(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


BLAST_DIR = first_existing_dir(
    ROOT / "resources" / "suppl20" / "blast_samples",
    ROOT / "Suppl20" / "blast_samples",
)
BLAST_RESULTS_DIR = first_existing_dir(
    ROOT / "resources" / "suppl20" / "blast_results",
    ROOT / "Suppl20" / "blast_results",
)
OUT_DIR = Path(__file__).parent / "suppl20a_plots"
OUT_DIR.mkdir(exist_ok=True)
SUMMARY_CSV = OUT_DIR / "suppl20a_quality_comparison_summary.csv"


def count_sequences_with_hits(blast_output_file: Path) -> tuple[int, int]:
    if not blast_output_file.exists():
        return 0, 0
    unique_queries_with_hits = set()
    total_hits = 0
    for line in blast_output_file.read_text().splitlines():
        if line.strip():
            query_id = line.split("\t")[0]
            unique_queries_with_hits.add(query_id)
            total_hits += 1
    return len(unique_queries_with_hits), total_hits


def analyze_homo_sapiens_hits(blast_output_file: Path) -> dict:
    if not blast_output_file.exists():
        return {
            "total_hs_queries": set(),
            "hs_mrna_queries": set(),
            "hs_rrna_queries": set(),
            "hs_other_queries": set(),
        }

    query_hits: dict[str, list[str]] = {}
    for line in blast_output_file.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        query_id = parts[0]
        subject_title = parts[12].lower()
        if "homo sapiens" not in subject_title:
            continue
        hit_type = "other"
        if (
            "ribosomal rna" in subject_title
            or "rrna" in subject_title
            or "18s" in subject_title
            or "28s" in subject_title
            or "5.8s" in subject_title
            or "5s ribosomal" in subject_title
        ):
            hit_type = "rrna"
        elif (
            "mrna" in subject_title
            or "messenger rna" in subject_title
            or "transcript variant" in subject_title
        ):
            hit_type = "mrna"
        query_hits.setdefault(query_id, []).append(hit_type)

    results = {
        "total_hs_queries": set(query_hits),
        "hs_mrna_queries": set(),
        "hs_rrna_queries": set(),
        "hs_other_queries": set(),
    }
    for query_id, hits in query_hits.items():
        if "rrna" in hits:
            results["hs_rrna_queries"].add(query_id)
        elif "mrna" in hits:
            results["hs_mrna_queries"].add(query_id)
        else:
            results["hs_other_queries"].add(query_id)
    return results


def fasta_query_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for line in path.read_text().splitlines() if line.startswith(">"))


def load_samples() -> list[str]:
    cache_file = BLAST_RESULTS_DIR / "blast_cache.json"
    if cache_file.exists():
        cache = json.loads(cache_file.read_text())
        samples = sorted({key.rsplit("_", 1)[0] for key in cache if key.endswith(("accepted", "rejected"))})
        if samples:
            return samples
    fasta_samples = sorted({path.name.replace("_accepted_overrepresented.fasta", "").replace("_rejected_overrepresented.fasta", "") for path in BLAST_DIR.glob("*_overrepresented.fasta") if "_accepted_" in path.name or "_rejected_" in path.name})
    return fasta_samples


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "svg"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


samples = load_samples()
print(f"Samples found: {samples}")

if not samples:
    raise FileNotFoundError(
        "No Suppl20 BLAST comparison samples were found. "
        f"Checked {BLAST_DIR} and {BLAST_RESULTS_DIR}."
    )

rows = []
canonical_sample = None
for sample in samples:
    sample_data = {
        "accepted": {"total": 0, "with_hits": 0, "no_hits": 0, "hs_mrna": 0, "hs_rrna": 0, "hs_other": 0, "non_hs": 0},
        "rejected": {"total": 0, "with_hits": 0, "no_hits": 0, "hs_mrna": 0, "hs_rrna": 0, "hs_other": 0, "non_hs": 0},
    }
    for label in ("accepted", "rejected"):
        fasta_file = BLAST_DIR / f"{sample}_{label}_overrepresented.fasta"
        blast_txt = BLAST_RESULTS_DIR / f"{sample}_{label}_blast.txt"
        total = fasta_query_count(fasta_file)
        with_hits, _ = count_sequences_with_hits(blast_txt)
        hs_results = analyze_homo_sapiens_hits(blast_txt)
        sample_data[label]["total"] = total
        sample_data[label]["with_hits"] = with_hits
        sample_data[label]["no_hits"] = total - with_hits
        sample_data[label]["hs_mrna"] = len(hs_results["hs_mrna_queries"])
        sample_data[label]["hs_rrna"] = len(hs_results["hs_rrna_queries"])
        sample_data[label]["hs_other"] = len(hs_results["hs_other_queries"])
        sample_data[label]["non_hs"] = with_hits - len(hs_results["total_hs_queries"])
        rows.append(
            {
                "sample": sample,
                "label": label,
                **sample_data[label],
            }
        )

    if sample_data["accepted"]["total"] == 0 and sample_data["rejected"]["total"] == 0:
        continue

    if canonical_sample is None:
        canonical_sample = sample

    categories = ["Accepted", "Rejected"]
    x = np.arange(len(categories))
    width = 0.5
    acc_total = sample_data["accepted"]["total"]
    rej_total = sample_data["rejected"]["total"]
    acc_nohits_pct = (sample_data["accepted"]["no_hits"] / acc_total * 100) if acc_total else 0
    rej_nohits_pct = (sample_data["rejected"]["no_hits"] / rej_total * 100) if rej_total else 0

    fig, ax = plt.subplots(figsize=(3, 3))
    fig.suptitle(f"Relative Composition: {sample}", fontsize=10, fontweight="bold")

    acc_mrna_pct = (sample_data["accepted"]["hs_mrna"] / acc_total * 100) if acc_total else 0
    acc_rrna_pct = (sample_data["accepted"]["hs_rrna"] / acc_total * 100) if acc_total else 0
    acc_hs_other_pct = (sample_data["accepted"]["hs_other"] / acc_total * 100) if acc_total else 0
    acc_non_hs_pct = (sample_data["accepted"]["non_hs"] / acc_total * 100) if acc_total else 0
    rej_mrna_pct = (sample_data["rejected"]["hs_mrna"] / rej_total * 100) if rej_total else 0
    rej_rrna_pct = (sample_data["rejected"]["hs_rrna"] / rej_total * 100) if rej_total else 0
    rej_hs_other_pct = (sample_data["rejected"]["hs_other"] / rej_total * 100) if rej_total else 0
    rej_non_hs_pct = (sample_data["rejected"]["non_hs"] / rej_total * 100) if rej_total else 0

    mrna_pct = [acc_mrna_pct, rej_mrna_pct]
    rrna_pct = [acc_rrna_pct, rej_rrna_pct]
    hs_other_pct = [acc_hs_other_pct, rej_hs_other_pct]
    non_hs_pct = [acc_non_hs_pct, rej_non_hs_pct]
    nohit_pct = [acc_nohits_pct, rej_nohits_pct]
    ax.bar(x, mrna_pct, width, label="H. sapiens mRNA", color="#2ecc71")
    ax.bar(x, rrna_pct, width, bottom=mrna_pct, label="H. sapiens rRNA", color="#9b59b6")
    bottom2 = np.array(mrna_pct) + np.array(rrna_pct)
    ax.bar(x, hs_other_pct, width, bottom=bottom2, label="H. sapiens other", color="#f39c12")
    bottom3 = bottom2 + np.array(hs_other_pct)
    ax.bar(x, non_hs_pct, width, bottom=bottom3, label="Non-H. sapiens", color="#95a5a6")
    bottom4 = bottom3 + np.array(non_hs_pct)
    ax.bar(x, nohit_pct, width, bottom=bottom4, label="No BLAST hits", color="#e74c3c")
    ax.set_ylabel("Percentage (%)", fontsize=8)
    ax.set_title("Relative Composition", fontsize=8, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_ylim(0, 100)
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(axis="y", alpha=0.3)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    save(fig, f"{sample}_quality_comparison")
    if sample == canonical_sample:
        save(fig, "suppl20a_quality_comparison")
    plt.close(fig)

summary_df = pd.DataFrame(rows)
summary_df.to_csv(SUMMARY_CSV, index=False)
print(f"Saved: {SUMMARY_CSV}")
print("Done.")
