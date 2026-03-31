"""
Figure 2a / Supplementary Figure 11b — Benchmarking notebook outputs
Faithful standalone reproduction of Suppl11/benchmarking_simplified_final.ipynb.
Figure 2a corresponds to the RNA amount bar plot; Supplementary Figure 11b
corresponds to the averaged SN vs Lysate hexbin figure.
Outputs to: refactoring_roadmap/suppl11b_2a_plots/

Notebook helper outputs remain available in code, but are disabled by default:
  2a_total_rna_amount_barplot.{svg,png}
  suppl11b_averaged_hexbin.{svg,png}
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir

ROOT = Path(__file__).parent.parent
DATA_DIR = get_generated_figure2_loaded_dir()
OUT_DIR = Path(__file__).parent / "suppl11b_2a_plots"
OUT_DIR.mkdir(exist_ok=True)

PLOT_FONT_SIZE = 8
PLOT_WIDTH = 3.5
PLOT_HEIGHT = 3.5
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = PLOT_FONT_SIZE

NAME_MAPPING = {
    "S1": "NTVE",
    "S2": "HIV-gag WT",
    "S3": "TRACE-seq",
    "S4": "MMLV-gag WT",
    "S5": "HIV-gag:MCP",
    "S6": "NLuc",
}
CONSTRUCTS = ["S1", "S2", "S3", "S4", "S5"]


def get_replicates(df: pd.DataFrame, construct: str, sample_type: str) -> pd.DataFrame:
    cols = [col for col in df.columns if col.startswith(f"{sample_type}_{construct}_")]
    return df[cols]


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


def plot_hexbin(ax: plt.Axes, x_data: pd.Series, y_data: pd.Series, title: str, is_mt: pd.Series | None = None):
    x = np.log10(x_data + 1)
    y = np.log10(y_data + 1)
    mask = (x > 0) & (y > 0)
    hb = ax.hexbin(x, y, gridsize=50, cmap="inferno", mincnt=1, bins="log", alpha=0.8, rasterized=True)
    if is_mt is not None:
        mt_mask = is_mt.values.astype(bool)
        if mt_mask.sum() > 0:
            ax.scatter(
                x[mt_mask], y[mt_mask], color="red", s=80, alpha=0.8,
                edgecolors="darkred", linewidths=1.5, marker="o", zorder=12, label="MT genes"
            )
    if mask.sum() > 10:
        corr, pval = pearsonr(x[mask], y[mask])
        ax.text(
            0.05, 0.95, f"r={corr:.3f}\np={pval:.2e}", transform=ax.transAxes,
            verticalalignment="top", fontsize=PLOT_FONT_SIZE + 1,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Lysate log10(TPM+1)", fontweight="bold")
    ax.set_ylabel("SN log10(TPM+1)", fontweight="bold")
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5, zorder=0, linewidth=1.5)
    return hb


print("Loading benchmarking data...")
tpm_df = pd.read_csv(DATA_DIR / "elo_tpm_matrix_protein_coding.csv", index_col=0)
gene_names = tpm_df["GeneName"]
tpm_df = tpm_df.drop(columns=["GeneName"])

averaged_data = pd.DataFrame(index=tpm_df.index)
averaged_data["GeneName"] = gene_names
for construct in CONSTRUCTS:
    sn_reps = get_replicates(tpm_df, construct, "SN")
    lys_reps = get_replicates(tpm_df, construct, "Lysate")
    averaged_data[f"SN_{construct}_avg"] = sn_reps.mean(axis=1)
    averaged_data[f"Lysate_{construct}_avg"] = lys_reps.mean(axis=1)

from ntvetools import load_gtf_df
gtf_data = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_is_mt = gtf_data["gene_id_is_mt"]
averaged_data["is_mt"] = [gene_id_is_mt.get(gene_id) for gene_id in averaged_data.index]

# Notebook helper outputs intentionally disabled by default:
#
# fig = plt.figure(figsize=(14, 2.5))
# gs = GridSpec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.15], wspace=0.3)
# axes = [fig.add_subplot(gs[0, i]) for i in range(5)]
# cbar_ax = fig.add_subplot(gs[0, 5])
#
# last_hb = None
# for i, construct in enumerate(CONSTRUCTS):
#     sn_col = f"SN_{construct}_avg"
#     lys_col = f"Lysate_{construct}_avg"
#     title = NAME_MAPPING[construct]
#     last_hb = plot_hexbin(axes[i], averaged_data[lys_col], averaged_data[sn_col], title, is_mt=averaged_data["is_mt"])
#
# fig.colorbar(last_hb, cax=cbar_ax, label="Count (log)")
# handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="red", markeredgecolor="darkred", markersize=8, label="MT genes")]
# axes[0].legend(handles=handles, loc="lower right", fontsize=PLOT_FONT_SIZE)
# save(fig, "suppl11b_averaged_hexbin")
# plt.close(fig)

numreads_df = pd.read_csv(DATA_DIR / "elo_numreads_matrix_protein_coding.csv", index_col=0)
mgl_reads = numreads_df[averaged_data["GeneName"] == "mGreenLantern"]
total_reads = numreads_df.sum()
mgl_fraction = mgl_reads.iloc[0] / total_reads
estimated_input_ng = (1.0 / mgl_fraction) / 1000.0
sn_cols = [col for col in estimated_input_ng.index if col.startswith("SN_")]
sn_input = estimated_input_ng[sn_cols]

mgl_plot_data = []
for construct in CONSTRUCTS:
    reps = [col for col in sn_cols if f"_{construct}_" in col]
    vals = sn_input[reps]
    for value in vals:
        mgl_plot_data.append({"Condition": NAME_MAPPING[construct], "Input_ng": value})

mgl_plot_df = pd.DataFrame(mgl_plot_data)
stats_df = mgl_plot_df.groupby("Condition")["Input_ng"].agg(["mean", "std"]).reset_index()
stats_df["Condition"] = pd.Categorical(stats_df["Condition"], categories=[NAME_MAPPING[c] for c in CONSTRUCTS], ordered=True)
stats_df = stats_df.sort_values("Condition")

# fig, ax = plt.subplots(figsize=(PLOT_WIDTH * 0.75, PLOT_HEIGHT))
# x_pos = np.arange(len(stats_df))
# bar_width = 0.55
# sn_color = "#3498db"
# ax.bar(x_pos, stats_df["mean"], bar_width, yerr=stats_df["std"], capsize=5, color=sn_color, alpha=0.7, error_kw={"linewidth": 1.5})
# for i, condition in enumerate(stats_df["Condition"]):
#     pts = mgl_plot_df[mgl_plot_df["Condition"] == condition]["Input_ng"]
#     jitter = np.random.normal(0, 0.05, size=len(pts))
#     ax.scatter(np.ones(len(pts)) * i + jitter, pts, color="black", s=15, alpha=0.6, zorder=10)
# ax.set_yscale("log")
# ax.axhline(1, ls="--", color="gray", alpha=0.5, linewidth=1.5, label="1 ng reference")
# ax.set_title(
#     "Benchmarking Experiment - Estimated RNA Input (SN only)\n(Bars: Mean ± SD, Points: Individual Replicates)",
#     fontweight="bold", fontsize=PLOT_FONT_SIZE + 1,
# )
# ax.set_ylabel("Estimated RNA Input (ng)", fontweight="bold")
# ax.set_xlabel("Condition", fontweight="bold")
# ax.set_xticks(x_pos)
# ax.set_xticklabels(stats_df["Condition"], rotation=45, ha="right")
# ax.legend(fontsize=PLOT_FONT_SIZE)
# ax.grid(True, alpha=0.3, axis="y")
# fig.tight_layout()
# save(fig, "2a_total_rna_amount_barplot")
# plt.close(fig)
print("Done.")
