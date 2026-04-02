"""
Supplementary Figure 12a — Benchmarking comparison CV of export ratio
Notebook-faithful standalone reproduction of the Elo combined overlay from
Figure2/cv_export_ratio.ipynb, with professional Benchmarking naming.
Outputs to: refactoring_roadmap/suppl12a_plots/
"""

from collections import defaultdict
from pathlib import Path
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys

FIG2 = Path(__file__).parent.parent / "Figure2"
sys.path.insert(0, str(FIG2))
sys.path.insert(0, str(Path(__file__).parent.parent))

from ntvetools import load_gtf_df
from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir
from rediscovery_analysis_lib import reorganize_samples_by_number

warnings.filterwarnings("ignore")
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 8

ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent / "suppl12a_plots"
OUT_DIR.mkdir(exist_ok=True)

plot_dir = OUT_DIR

name_mapping = {
    "Benchmarking-S1": "NTVE",
    "Benchmarking-S2": "HIV-gag WT",
    "Benchmarking-S3": "TRACE-seq",
    "Benchmarking-S4": "MMLV-gag WT",
    "Benchmarking-S5": "HIV-gag:MCP",
    "Benchmarking-S6": "NLuc",
}


def get_sample_display_name(exp_prefix, sample_num):
    key = f"{exp_prefix}-S{sample_num}"
    mapped_name = name_mapping.get(key)
    if mapped_name:
        return f"S{sample_num}: {mapped_name}"
    return f"S{sample_num}"


def save_svg(fig, filename):
    output_path = plot_dir / filename
    fig.savefig(output_path, format="svg", bbox_inches="tight")
    png_path = output_path.with_suffix(".png")
    fig.savefig(png_path, format="png", bbox_inches="tight", dpi=150)
    print(f"Saved: {output_path}")
    print(f"Saved: {png_path}")


def aggregate_transcripts_to_genes(tpm_df, transcript_info_df):
    transcript_to_gene = transcript_info_df.set_index("transcript_id")["gene_id"].to_dict()
    tpm_with_gene = tpm_df.copy()
    tpm_with_gene["gene_id"] = tpm_with_gene.index.map(lambda t: transcript_to_gene.get(t, ""))
    return tpm_with_gene.groupby("gene_id").sum(numeric_only=True)


def create_gene_df(gene_tpm, gene_numreads, gtf_df_data, gene_id_to_gene_name, gene_id_is_mt):
    gtf_df = gtf_df_data["gtf_df"]
    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype = gene_level_gtf.set_index("gene_id")
    gene_biotype_dict = defaultdict(lambda: "NA", gene_biotype["gene_biotype"].to_dict())

    gene_df = gene_tpm.copy()
    gene_df.index.name = "GeneID"
    gene_df["GeneName"] = gene_df.index.map(lambda g: gene_id_to_gene_name.get(g, g))
    gene_df["GeneBiotype"] = gene_df.index.map(lambda g: gene_biotype_dict.get(g, "NA"))
    gene_df["is_MT"] = gene_df.index.map(lambda g: gene_id_is_mt.get(g, False))

    for col in gene_numreads.columns:
        gene_df[col.replace("_NumReads", "") + "_NumReads"] = gene_numreads[col]

    return gene_df


def filter_genes(gene_df, min_biotype="protein_coding", exclude_mt=True, exclude_custom=True):
    df = gene_df.copy()

    if min_biotype:
        df = df[df["GeneBiotype"] == min_biotype].copy()
        print(f"After {min_biotype} filter: {df.shape[0]} genes")

    if exclude_mt:
        df = df[df["is_MT"] == False].copy()
        print(f"After removing MT genes: {df.shape[0]} genes")

    if exclude_custom:
        custom_genes = ["mGreenLantern", "Gag_common_core", "HIV_gag_MCP", "MMLV_gag", "minigag_PABP", "TetOn3G-V9I_mScarlet", "eUnaG2"]
        df = df[~df["GeneName"].isin(custom_genes)].copy()
        print(f"After removing custom genes: {df.shape[0]} genes")

    return df


def renormalize_tpm(df):
    df_copy = df.copy()
    numeric_cols = [col for col in df_copy.columns if col not in ["GeneName", "GeneBiotype", "is_MT"]]
    numeric_cols = [col for col in numeric_cols if not col.endswith("_NumReads")]

    for col in numeric_cols:
        tpm_sum = df_copy[col].sum()
        if tpm_sum > 0:
            df_copy[col] = df_copy[col] * 1e6 / tpm_sum

    return df_copy


def parse_sample_names(numeric_cols_elo):
    sample_name_mapping = {}
    for col in numeric_cols_elo:
        match = re.match(r"SN_S(\d+)_(\d+)", str(col))
        if match:
            sample_name_mapping[col] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "SN"}
            continue
        match = re.match(r"Lysate_S(\d+)_(\d+)", str(col))
        if match:
            sample_name_mapping[col] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "Lysate"}
    return sample_name_mapping


def compute_export_ratio_cv(organized_sn, organized_lysate, gene_df, sample_name_mapping):
    results = {}

    for sample_num in sorted(organized_sn.keys()):
        if sample_num not in organized_lysate:
            continue

        sn_reps = organized_sn[sample_num]
        lys_reps = organized_lysate[sample_num]

        sn_by_rep = {}
        for col in sn_reps:
            rep_num = sample_name_mapping[col]["Replicate"]
            sn_by_rep[rep_num] = col

        lys_by_rep = {}
        for col in lys_reps:
            rep_num = sample_name_mapping[col]["Replicate"]
            lys_by_rep[rep_num] = col

        common_reps = sorted(set(sn_by_rep.keys()) & set(lys_by_rep.keys()))
        if len(common_reps) < 2:
            print(f"  Sample {sample_num}: only {len(common_reps)} matched replicates, skipping")
            continue

        lysate_cols = [lys_by_rep[rep_num] for rep_num in common_reps]
        valid_mask = (gene_df[lysate_cols] > 0).all(axis=1)

        ratio_df = pd.DataFrame(index=gene_df.index)
        for rep_num in common_reps:
            sn_col = sn_by_rep[rep_num]
            lys_col = lys_by_rep[rep_num]
            ratio_df.loc[valid_mask, f"ratio_rep{rep_num}"] = (
                gene_df.loc[valid_mask, sn_col] / gene_df.loc[valid_mask, lys_col]
            )

        ratio_cols = [c for c in ratio_df.columns if c.startswith("ratio_")]

        result = pd.DataFrame(index=gene_df.index)
        result["mean_lysate_tpm"] = gene_df[lysate_cols].mean(axis=1)
        result["mean_ratio"] = ratio_df[ratio_cols].mean(axis=1)
        result["std_ratio"] = ratio_df[ratio_cols].std(axis=1)
        result["cv_ratio"] = result["std_ratio"] / result["mean_ratio"]
        result["cv_ratio"] = result["cv_ratio"].replace([np.inf, -np.inf], np.nan)
        result["n_reps"] = len(common_reps)

        results[sample_num] = result
        print(
            f"  Sample {sample_num}: {len(common_reps)} matched replicates, "
            f"{int(valid_mask.sum())} genes with nonzero lysate in all matched replicates, "
            f"{result['cv_ratio'].notna().sum()} genes with valid CV"
        )

    return results


def compute_cv_quantile_curve(cv_result, direction="bottom", n_quantiles=100):
    df = cv_result.dropna(subset=["cv_ratio"]).copy()
    df = df[df["mean_lysate_tpm"] > 1.0]
    df = df.sort_values("mean_lysate_tpm")

    n_genes = len(df)
    quantiles = np.linspace(1, 100, n_quantiles)
    median_cvs = []

    for q in quantiles:
        n_select = max(1, int(np.ceil(n_genes * q / 100)))
        if direction == "bottom":
            subset = df.iloc[:n_select]
        else:
            subset = df.iloc[-n_select:]
        median_cvs.append(subset["cv_ratio"].median())

    return quantiles, np.array(median_cvs), n_genes


print("Loading GTF reference data...")
gtf_df_data = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_gene_name_def = defaultdict(str, gtf_df_data["gene_id_to_name"])
gene_id_is_mt = gtf_df_data.get("gene_id_is_mt", {})

output_dir = get_generated_figure2_loaded_dir()
print("Loading gene-level data...")
elo_tpm = pd.read_csv(output_dir / "elo_transcript_tpm_matrix.csv", index_col=0)
elo_numreads = pd.read_csv(output_dir / "elo_transcript_numreads_matrix.csv", index_col=0)
elo_transcript_info = pd.read_csv(output_dir / "elo_transcript_info.csv")
print(f"Elo: {elo_tpm.shape[0]} transcripts, {elo_tpm.shape[1]} samples")

elo_gene_tpm = aggregate_transcripts_to_genes(elo_tpm, elo_transcript_info)
elo_gene_numreads = aggregate_transcripts_to_genes(elo_numreads, elo_transcript_info)
print(f"Elo gene-level: {elo_gene_tpm.shape[0]} genes")

gene_df_elo = create_gene_df(elo_gene_tpm, elo_gene_numreads, gtf_df_data, gene_id_to_gene_name_def, gene_id_is_mt)
print(f"Elo: {gene_df_elo.shape[0]} genes (before filtering)")

print("\nFiltering Benchmarking:")
gene_df_elo_pc = filter_genes(gene_df_elo)

print("Renormalizing TPM...")
gene_df_elo_pc = renormalize_tpm(gene_df_elo_pc)
numeric_cols_elo = [col for col in gene_df_elo_pc.columns if col not in ["GeneName", "GeneBiotype", "is_MT"] and not col.endswith("_NumReads")]
print(f"Elo samples: {len(numeric_cols_elo)}")

sample_name_mapping = parse_sample_names(numeric_cols_elo)
experiments = {
    "Benchmarking": {
        "SN": [c for c, meta in sample_name_mapping.items() if meta["Type"] == "SN"],
        "Lysate": [c for c, meta in sample_name_mapping.items() if meta["Type"] == "Lysate"],
    }
}

elo_sn_samples = experiments["Benchmarking"]["SN"]
elo_lysate_samples = experiments["Benchmarking"]["Lysate"]
elo_organized = {
    "SN": reorganize_samples_by_number(elo_sn_samples, sample_name_mapping),
    "Lysate": reorganize_samples_by_number(elo_lysate_samples, sample_name_mapping),
}

print("Benchmarking sample organization:")
for sample_type in elo_organized:
    print(f"\n  {sample_type}:")
    for sample_num in sorted(elo_organized[sample_type].keys()):
        reps = elo_organized[sample_type][sample_num]
        print(f"    Sample {sample_num}: {len(reps)} replicates")

elo_colors = {
    1: "#77C8DC",
    2: "#BDC7E3",
    3: "#ECB0A3",
    4: "#D7E6D5",
    5: "#B4C7B0",
    6: "#9E9E9E",
}
elo_linestyles = {
    1: "-",
    2: (0, (8, 3)),
    3: (0, (3, 2)),
    4: "-.",
    5: (0, (8, 3, 1, 3, 1, 3)),
    6: (0, (8, 3)),
}
elo_linewidths = {
    1: 3.5,
    2: 2.0,
    3: 2.0,
    4: 2.0,
    5: 2.0,
    6: 2.0,
}

print("Computing export ratio CV for Benchmarking experiment...")
elo_cv_results = compute_export_ratio_cv(
    elo_organized["SN"], elo_organized["Lysate"], gene_df_elo_pc, sample_name_mapping
)

elo_cv_curves = {}
for sample_num in sorted(elo_cv_results.keys()):
    elo_cv_curves[sample_num] = {}
    for direction in ["bottom", "top"]:
        quantiles, median_cvs, n_genes = compute_cv_quantile_curve(elo_cv_results[sample_num], direction=direction)
        elo_cv_curves[sample_num][direction] = {
            "quantiles": quantiles,
            "median_cvs": median_cvs,
            "n_genes": n_genes,
        }
    print(
        f"  Sample {sample_num}: {n_genes} genes, "
        f"median CV (bottom 50%)={elo_cv_curves[sample_num]['bottom']['median_cvs'][49]:.3f}, "
        f"median CV (top 50%)={elo_cv_curves[sample_num]['top']['median_cvs'][49]:.3f}"
    )

# Combined overlay: 1x2 panels matching the manuscript legend
fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
panel_meta = [
    ("bottom", "Lowest to highest expression"),
    ("top", "Highest to lowest expression"),
]

for dir_idx, (direction, panel_title) in enumerate(panel_meta):
    ax = axes[dir_idx]
    for sample_num in sorted(elo_cv_curves.keys()):
        curve = elo_cv_curves[sample_num][direction]
        cv_at_50 = curve["median_cvs"][49]
        label = f"{get_sample_display_name('Benchmarking', sample_num)} (CV@50%={cv_at_50:.2f})"
        ax.plot(
            curve["quantiles"],
            curve["median_cvs"],
            label=label,
            color=elo_colors[sample_num],
            linewidth=elo_linewidths[sample_num],
            linestyle=elo_linestyles[sample_num],
            alpha=0.95,
        )

    ax.set_xlim(-5, 105)
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("Cumulative percentile of detected genes (%)", fontweight="bold", fontsize=8)
    if dir_idx == 0:
        ax.set_ylabel("Coefficient of variation of export ratio", fontweight="bold", fontsize=8)
    ax.set_title(panel_title, fontweight="bold", fontsize=9)
    ax.legend(loc="best", fontsize=7, frameon=True, framealpha=0.85)

fig.suptitle("Benchmarking comparison - coefficient of variation across cumulative expression percentiles", fontsize=10, fontweight="bold", y=1.02)
plt.tight_layout()
save_svg(fig, "suppl12a_benchmarking_cv_export_ratio_combined.svg")
plt.close(fig)
print("Done.")
