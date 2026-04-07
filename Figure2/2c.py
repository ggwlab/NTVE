"""
Figure 2c — Benchmarking comparison transcript-length correlation
Notebook-faithful standalone reproduction of the matched, non-bootstrap Elo panel
from Figure2/correlation_by_transcript_length_noMT_1000bp.ipynb.
Outputs to: refactoring_roadmap/2c_plots/
"""

import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df
from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir

warnings.filterwarnings("ignore")

# Global plot parameters
PLOT_FONT_SIZE = 8
PLOT_WIDTH = 12
PLOT_HEIGHT = 2.5
PLOT_OUTPUT_DIR = str(Path(__file__).parent / "2c_plots")
CSV_OUTPUT_DIR = Path(__file__).parent / "2c_csv"

# Set matplotlib for SVG editability - no embedded fonts
plt.rcParams["svg.fonttype"] = "none"  # Don't convert text to paths
plt.rcParams["font.size"] = PLOT_FONT_SIZE

# Create output directory if it doesn't exist
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)
CSV_OUTPUT_DIR.mkdir(exist_ok=True)

print("Imports successful!")
print(f"✓ Plot output directory: {PLOT_OUTPUT_DIR}")


def create_biotype_dicts(transcript_info_df, gtf_df_data):
    """Create mapping from transcript_id to biotype (protein_coding)"""
    gtf_df = gtf_df_data["gtf_df"]
    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype_dict = defaultdict(lambda: "NA", gene_level_gtf.set_index("gene_id")["gene_biotype"].to_dict())

    transcript_to_gene = transcript_info_df.set_index("transcript_id")["gene_id"].to_dict()
    transcript_biotype_dict = {}
    for transcript_id, gene_id in transcript_to_gene.items():
        transcript_biotype_dict[transcript_id] = gene_biotype_dict.get(gene_id, "NA")

    return transcript_biotype_dict


def parse_sample_names_transcript(tpm_df, lims_mapping=None):
    """Parse sample names and create mapping with transcript length info"""
    sample_name_mapping = {}

    for col_name in tpm_df.columns:
        col_str = str(col_name)

        if lims_mapping is not None and col_str in lims_mapping:
            lims_name = lims_mapping[col_str]
            match = re.match(r"SN_(\d+)_(\d+)", str(lims_name))
            if match:
                sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "SN", "Experiment": "Spam", "ColumnName": col_str}
                continue
            match = re.match(r"L_(\d+)_(\d+)", str(lims_name))
            if match:
                sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "L", "Experiment": "Spam", "ColumnName": col_str}

        match = re.match(r"SN_S(\d+)_(\d+)", col_str)
        if match:
            sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "SN", "Experiment": "Benchmarking", "ColumnName": col_str}
            continue
        match = re.match(r"Lysate_S(\d+)_(\d+)", col_str)
        if match:
            sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "Lysate", "Experiment": "Benchmarking", "ColumnName": col_str}

    return sample_name_mapping


def organize_samples_by_type_and_number(sample_map, tpm_df):
    """Organize samples by type and sample number"""
    organized = {}
    for col_name, info in sample_map.items():
        if col_name not in tpm_df.columns:
            continue
        sample_type = info["Type"]
        sample_num = info["Sample"]
        if sample_type not in organized:
            organized[sample_type] = {}
        if sample_num not in organized[sample_type]:
            organized[sample_type][sample_num] = []
        organized[sample_type][sample_num].append(col_name)
    return organized


def create_length_bins():
    """Create transcript length bins: 0-999, 1000-1999, ..., 4000-4999, 5000+"""
    bins = [0, 1000, 2000, 3000, 4000, 5000, float("inf")]
    labels = ["0-999", "1000-1999", "2000-2999", "3000-3999", "4000-4999", "5000+"]
    return bins, labels


def get_sample_pairs(source_samples, target_samples, sample_map, pairing_mode="all_pairs"):
    if pairing_mode == "all_pairs":
        return [(source_col, target_col) for source_col in source_samples for target_col in target_samples]

    if pairing_mode == "matched_replicates":
        source_by_replicate = defaultdict(list)
        target_by_replicate = defaultdict(list)

        for source_col in source_samples:
            replicate = sample_map.get(source_col, {}).get("Replicate")
            if replicate is not None:
                source_by_replicate[replicate].append(source_col)

        for target_col in target_samples:
            replicate = sample_map.get(target_col, {}).get("Replicate")
            if replicate is not None:
                target_by_replicate[replicate].append(target_col)

        sample_pairs = []
        for replicate in sorted(set(source_by_replicate.keys()) & set(target_by_replicate.keys())):
            for source_col in source_by_replicate[replicate]:
                for target_col in target_by_replicate[replicate]:
                    sample_pairs.append((source_col, target_col))
        return sample_pairs

    raise ValueError(f"Unknown pairing_mode: {pairing_mode}")


def analyze_correlation_by_length(tpm_df, sample_map, transcript_info_df, biotype_dict, organized, source_type, target_type, tpm_threshold=1.0, min_biotype="protein_coding", pairing_mode="all_pairs"):
    bins, bin_labels = create_length_bins()
    pc_transcripts = [tid for tid, biotype in biotype_dict.items() if biotype == min_biotype and tid in tpm_df.index]
    tpm_pc = tpm_df.loc[pc_transcripts].copy()
    length_dict = transcript_info_df.set_index("transcript_id")["length"].to_dict()
    tpm_pc["length"] = tpm_pc.index.map(lambda t: length_dict.get(t, np.nan))
    tpm_pc = tpm_pc.dropna(subset=["length"])
    tpm_pc["length_bin"] = pd.cut(tpm_pc["length"], bins=bins, labels=bin_labels, right=False)

    results = {}
    sample_numbers = sorted(organized.get(source_type, {}).keys())

    for sample_num in sample_numbers:
        source_samples = organized[source_type].get(sample_num, [])
        target_samples = organized[target_type].get(sample_num, [])
        results[sample_num] = {}

        # Filter samples to only those present in the data
        source_samples = [s for s in source_samples if s in tpm_pc.columns]
        target_samples = [t for t in target_samples if t in tpm_pc.columns]

        sample_pairs = get_sample_pairs(source_samples, target_samples, sample_map, pairing_mode=pairing_mode)

        if not source_samples or not target_samples or not sample_pairs:
            for bin_label in bin_labels:
                results[sample_num][bin_label] = {"pearson_r_mean": np.nan, "pearson_r_std": np.nan, "pearson_pval_mean": np.nan, "spearman_r_mean": np.nan, "spearman_r_std": np.nan, "spearman_pval_mean": np.nan, "n_transcripts": 0, "n_pairwise_correlations": 0, "pearson_r_values": [], "spearman_r_values": []}
            continue

        for bin_label in bin_labels:
            bin_data = tpm_pc[tpm_pc["length_bin"] == bin_label].copy()
            if len(bin_data) == 0:
                results[sample_num][bin_label] = {"pearson_r_mean": np.nan, "pearson_r_std": np.nan, "pearson_pval_mean": np.nan, "spearman_r_mean": np.nan, "spearman_r_std": np.nan, "spearman_pval_mean": np.nan, "n_transcripts": 0, "n_pairwise_correlations": 0, "pearson_r_values": [], "spearman_r_values": []}
                continue

            # Find transcripts expressed in ALL source AND ALL target replicates
            expressed_in_all_source = (bin_data[source_samples] > tpm_threshold).all(axis=1)
            expressed_in_all_target = (bin_data[target_samples] > tpm_threshold).all(axis=1)
            expressed_in_all = expressed_in_all_source & expressed_in_all_target

            n_expressed_all = expressed_in_all.sum()

            if n_expressed_all <= 1:
                results[sample_num][bin_label] = {"pearson_r_mean": np.nan, "pearson_r_std": np.nan, "pearson_pval_mean": np.nan, "spearman_r_mean": np.nan, "spearman_r_std": np.nan, "spearman_pval_mean": np.nan, "n_transcripts": n_expressed_all, "n_pairwise_correlations": 0, "pearson_r_values": [], "spearman_r_values": []}
                continue

            # Filter to only commonly expressed transcripts
            bin_data_filtered = bin_data.loc[expressed_in_all]

            pearson_rs, pearson_pvals, spearman_rs, spearman_pvals = [], [], [], []

            # Calculate correlations for each replicate pairing strategy
            for source_col, target_col in sample_pairs:
                source_values = bin_data_filtered[source_col].values
                target_values = bin_data_filtered[target_col].values

                pearson_r, pearson_pval = pearsonr(source_values, target_values)
                spearman_r, spearman_pval = spearmanr(source_values, target_values)
                pearson_rs.append(pearson_r)
                pearson_pvals.append(pearson_pval)
                spearman_rs.append(spearman_r)
                spearman_pvals.append(spearman_pval)

            n_pairwise = len(pearson_rs)
            pearson_std = np.std(pearson_rs) if n_pairwise > 2 else np.nan
            spearman_std = np.std(spearman_rs) if n_pairwise > 2 else np.nan
            results[sample_num][bin_label] = {"pearson_r_mean": np.mean(pearson_rs), "pearson_r_std": pearson_std, "pearson_pval_mean": np.mean(pearson_pvals), "spearman_r_mean": np.mean(spearman_rs), "spearman_r_std": spearman_std, "spearman_pval_mean": np.mean(spearman_pvals), "n_transcripts": int(n_expressed_all), "n_pairwise_correlations": n_pairwise, "pearson_r_values": [float(x) for x in pearson_rs], "spearman_r_values": [float(x) for x in spearman_rs]}

    return results, bin_labels


# Name mapping
name_mapping = {
    "Benchmarking-S1": "NTVE",
    "Benchmarking-S2": "HIV-gag WT",
    "Benchmarking-S3": "TRACE-seq",
    "Benchmarking-S4": "MMLV-gag WT",
    "Benchmarking-S5": "HIV-gag:MCP",
    "Benchmarking-S6": "NLuc",
}


def add_value_dots(ax, x_center, values, width, color):
    if not values:
        return
    if len(values) == 1:
        x_offsets = np.array([0.0])
    else:
        x_offsets = np.linspace(-width * 0.18, width * 0.18, len(values))
    ax.scatter(x_center + x_offsets, values, s=12, color=color, edgecolors="none", alpha=0.6, zorder=4)


def plot_correlation_by_length(results, bin_labels, experiment_name, source_type, target_type, save_filename=None):
    # 2×3 grid: row 0 = [S1, S2, S3], row 1 = [S5, S4, S6]
    layout = [[1, 2, 3], [5, 4, 6]]
    nrows, ncols = 2, 3

    panel_config = {
        1: dict(
            title="NTVE$_{PABP}$",
            banner="HIV-1 Gag\npoly(A)-tail / PABP",
            color="#7EB8D9",
            left_anno=True,
        ),
        2: dict(
            title="this publication",
            banner="HIV-1 Gag\nunspecific",
            color="#A39CC8",
            left_anno=False,
        ),
        3: dict(
            title="TRACE-seq",
            banner="endogenous EVs\n5\u2019 cap / EIF4E + YTHDF1",
            color="#E8927A",
            left_anno=False,
        ),
        5: dict(
            title="COURIER",
            banner="HIV-1 Gag-MCP\nunspecific/MCP",
            color="#74B874",
            left_anno=True,
        ),
        4: dict(
            title="COURIER",
            banner="MMLV Gag\nunspecific",
            color="#B3DBA3",
            left_anno=False,
        ),
        6: dict(
            title="transfection control",
            banner="transfection\ncontrol",
            color="#F4C48C",
            left_anno=False,
        ),
    }

    # First pass: determine global max for n_transcripts
    max_n_transcripts = 0
    for sample_num in results:
        for bin_label in bin_labels:
            if bin_label in results[sample_num]:
                n_transcripts = results[sample_num][bin_label].get("n_transcripts", 0)
                max_n_transcripts = max(max_n_transcripts, n_transcripts)

    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 9))
    plt.subplots_adjust(hspace=0.6, wspace=0.55, top=0.88, bottom=0.10, left=0.08, right=0.95)

    x_pos = np.arange(len(bin_labels))
    width = 0.6

    for row_idx, row in enumerate(layout):
        for col_idx, sample_num in enumerate(row):
            ax = axes[row_idx, col_idx]
            cfg = panel_config.get(sample_num, {})

            if sample_num not in results:
                ax.set_visible(False)
                continue

            spearman_means, spearman_stds, n_transcripts_list, spearman_value_sets = [], [], [], []
            for bin_label in bin_labels:
                if bin_label in results[sample_num]:
                    data = results[sample_num][bin_label]
                    spearman_means.append(data["spearman_r_mean"])
                    spearman_stds.append(data["spearman_r_std"])
                    n_transcripts_list.append(data.get("n_transcripts", 0))
                    spearman_value_sets.append(data.get("spearman_r_values", []))
                else:
                    spearman_means.append(np.nan)
                    spearman_stds.append(np.nan)
                    n_transcripts_list.append(0)
                    spearman_value_sets.append([])

            ax.bar(x_pos, spearman_means, width, alpha=0.8, color="darkgrey")

            for idx, values in enumerate(spearman_value_sets):
                if len(values) > 2 and np.isfinite(spearman_means[idx]) and np.isfinite(spearman_stds[idx]):
                    ax.errorbar(x_pos[idx], spearman_means[idx], yerr=spearman_stds[idx], fmt="none", ecolor="black", elinewidth=0.5, capsize=5, zorder=5)
                add_value_dots(ax, x_pos[idx], values, width, "black")

            ax.set_xlabel("Transcript Length", fontsize=PLOT_FONT_SIZE)
            ax.set_ylabel("Spearman correlation coefficient", fontsize=PLOT_FONT_SIZE)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(bin_labels, rotation=45, ha="right", fontsize=PLOT_FONT_SIZE - 1)
            ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax.tick_params(axis="y", labelsize=PLOT_FONT_SIZE - 1)
            ax.grid(True, alpha=0.3, axis="y")
            ax.set_ylim(0, 1.0)

            # Bold panel title inside (upper left)
            ax.text(0.05, 0.97, cfg.get("title", ""), transform=ax.transAxes,
                    fontsize=PLOT_FONT_SIZE + 1, fontweight="bold", va="top", ha="left")

            # Right y-axis: number of transcripts
            ax2 = ax.twinx()
            ax2.plot(x_pos, n_transcripts_list, color="darkgreen", marker="o",
                     linewidth=2, markersize=5, alpha=0.7)
            ax2.set_ylabel("Number of transcripts", fontsize=PLOT_FONT_SIZE, color="darkgreen")
            ax2.tick_params(axis="y", labelcolor="darkgreen", labelsize=PLOT_FONT_SIZE - 1)
            ax2.set_ylim(0, max_n_transcripts * 1.1)

            # Colored banner above the panel
            banner_text = cfg.get("banner", "")
            banner_color = cfg.get("color", "#DDDDDD")
            ax.annotate(
                banner_text,
                xy=(0.5, 1.0),
                xycoords="axes fraction",
                xytext=(0.5, 1.03),
                textcoords="axes fraction",
                ha="center", va="bottom",
                fontsize=PLOT_FONT_SIZE - 1,
                clip_on=False,
                annotation_clip=False,
                bbox=dict(boxstyle="round,pad=0.3", facecolor=banner_color, edgecolor="none", alpha=0.9),
            )

            # "budding module / RNA motif/binder" label for first column
            if cfg.get("left_anno", False):
                ax.annotate(
                    "budding module\nRNA motif/binder",
                    xy=(0.0, 1.0),
                    xycoords="axes fraction",
                    xytext=(-0.02, 1.10),
                    textcoords="axes fraction",
                    ha="right", va="bottom",
                    fontsize=PLOT_FONT_SIZE - 2,
                    clip_on=False,
                    annotation_clip=False,
                    color="#444444",
                )

    if save_filename:
        output_path = os.path.join(PLOT_OUTPUT_DIR, save_filename)
        fig.savefig(output_path, format="svg", bbox_inches="tight")
        png_path = os.path.splitext(output_path)[0] + ".png"
        fig.savefig(png_path, format="png", bbox_inches="tight", dpi=150)
        print(f"Saved: {output_path}")
        print(f"Saved: {png_path}")

    return fig


def export_correlation_tables(results, bin_labels, stem: str) -> None:
    summary_rows = []
    replicate_rows = []

    for sample_num in sorted(results.keys()):
        for bin_label in bin_labels:
            if bin_label not in results[sample_num]:
                continue
            data = results[sample_num][bin_label]
            summary_rows.append(
                {
                    "sample_num": sample_num,
                    "length_bin": bin_label,
                    "n_transcripts": data.get("n_transcripts", 0),
                    "n_pairwise_correlations": data.get("n_pairwise_correlations", 0),
                    "pearson_r_mean": data.get("pearson_r_mean"),
                    "pearson_r_std": data.get("pearson_r_std"),
                    "pearson_pval_mean": data.get("pearson_pval_mean"),
                    "spearman_r_mean": data.get("spearman_r_mean"),
                    "spearman_r_std": data.get("spearman_r_std"),
                    "spearman_pval_mean": data.get("spearman_pval_mean"),
                }
            )
            pearson_values = data.get("pearson_r_values", [])
            spearman_values = data.get("spearman_r_values", [])
            for pair_idx, (pearson_r, spearman_r) in enumerate(zip(pearson_values, spearman_values), start=1):
                replicate_rows.append(
                    {
                        "sample_num": sample_num,
                        "length_bin": bin_label,
                        "pair_index": pair_idx,
                        "pearson_r": pearson_r,
                        "spearman_r": spearman_r,
                    }
                )

    pd.DataFrame(summary_rows).to_csv(CSV_OUTPUT_DIR / f"{stem}_summary.csv", index=False)
    pd.DataFrame(replicate_rows).to_csv(CSV_OUTPUT_DIR / f"{stem}_pairwise_values.csv", index=False)


print("Loading transcript-level data...")
output_dir = get_generated_figure2_loaded_dir()
elo_tpm = pd.read_csv(output_dir / "elo_transcript_tpm_matrix.csv", index_col=0)
elo_transcript_info = pd.read_csv(output_dir / "elo_transcript_info.csv")
print(f"Elo: {elo_tpm.shape[0]} transcripts, {elo_tpm.shape[1]} samples")

print("Loading GTF reference data...")
gtf_df_data = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
elo_biotype_dict = create_biotype_dicts(elo_transcript_info, gtf_df_data)
print(f"Elo protein_coding transcripts: {sum(1 for b in elo_biotype_dict.values() if b == 'protein_coding')}")

elo_sample_map = parse_sample_names_transcript(elo_tpm)
print(f"Mapped {len(elo_sample_map)} Benchmarking samples")
elo_organized = organize_samples_by_type_and_number(elo_sample_map, elo_tpm)
print("Sample organization complete.")

bins, bin_labels = create_length_bins()
print(f"Created {len(bin_labels)} transcript length bins (1000bp each)")

print("\nAnalyzing Benchmarking experiment (matched replicates)...")
elo_corr_results_matched, elo_bins_matched = analyze_correlation_by_length(
    elo_tpm,
    elo_sample_map,
    elo_transcript_info,
    elo_biotype_dict,
    elo_organized,
    "SN",
    "Lysate",
    pairing_mode="matched_replicates",
)
print(f"✓ Processed {len(elo_corr_results_matched)} Benchmarking samples")

print("Creating Benchmarking correlation plots (1000bp bins, matched replicates)...")
export_correlation_tables(
    elo_corr_results_matched,
    elo_bins_matched,
    "2c_benchmarking_correlation_by_length_matched",
)
fig_elo_corr_matched = plot_correlation_by_length(
    elo_corr_results_matched,
    elo_bins_matched,
    "Benchmarking Experiment (Matched Replicates)",
    "SN",
    "Lysate",
    save_filename="2c_benchmarking_correlation_by_length_matched.svg",
)
plt.close(fig_elo_corr_matched)
print("Done.")
