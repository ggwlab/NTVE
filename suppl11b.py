"""
Supplementary Figure 11b — Benchmarking comparison averaged SN vs Lysate hexbins
Standalone reproduction of the hexbin panel from Suppl11/benchmarking_simplified_final.ipynb.
Outputs to: refactoring_roadmap/suppl11b_plots/
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from ntvetools import load_gtf_df
from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir

sns.set_style("whitegrid")

# Global plot parameters
PLOT_FONT_SIZE = 8
PLOT_WIDTH = 3.5
PLOT_HEIGHT = 3.5

# Set matplotlib for SVG editability - no embedded fonts
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = PLOT_FONT_SIZE

ROOT = Path(__file__).parent.parent
DATA_DIR = get_generated_figure2_loaded_dir()
OUT_DIR = Path(__file__).parent / "suppl11b_plots"
OUT_DIR.mkdir(exist_ok=True)

NAME_MAPPING = {
    "S1": "NTVE",
    "S2": "HIV-gag WT",
    "S3": "TRACE-seq",
    "S4": "MMLV-gag WT",
    "S5": "HIV-gag:MCP",
    "S6": "NLuc",
}

# Replicates to average
CONSTRUCTS = ["S1", "S2", "S3", "S4", "S5"]


def get_replicates(df: pd.DataFrame, construct: str, sample_type: str) -> pd.DataFrame:
    cols = [c for c in df.columns if c.startswith(f"{sample_type}_{construct}_")]
    return df[cols]


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


print("Loading benchmarking data...")
tpm_df = pd.read_csv(DATA_DIR / "elo_tpm_matrix_protein_coding.csv", index_col=0)
gene_names = tpm_df["GeneName"]
tpm_df = tpm_df.drop(columns=["GeneName"])


def plot_hexbin(ax: plt.Axes, x_data: pd.Series, y_data: pd.Series, title: str, is_mt: pd.Series | None = None):
    x = np.log10(x_data + 1)
    y = np.log10(y_data + 1)

    mask = (x > 0) & (y > 0)

    hb = ax.hexbin(x, y, gridsize=50, cmap="inferno", mincnt=1, bins="log", alpha=0.8, rasterized=True)

    # Overlay: Mitochondrial genes as red scatter points
    if is_mt is not None:
        mt_mask = is_mt.values.astype(bool)
        if mt_mask.sum() > 0:
            ax.scatter(
                x[mt_mask],
                y[mt_mask],
                color="red",
                s=80,
                alpha=0.8,
                edgecolors="darkred",
                linewidths=1.5,
                marker="o",
                zorder=12,
                label="MT genes",
            )

    if mask.sum() > 10:
        corr, pval = pearsonr(x[mask], y[mask])
        ax.text(
            0.05,
            0.95,
            f"r={corr:.3f}\np={pval:.2e}",
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=PLOT_FONT_SIZE + 1,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Lysate log10(TPM+1)", fontweight="bold")
    ax.set_ylabel("SN log10(TPM+1)", fontweight="bold")

    # Add diagonal
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5, zorder=0, linewidth=1.5)
    return hb


averaged_data = pd.DataFrame(index=tpm_df.index)
averaged_data["GeneName"] = gene_names
for construct in CONSTRUCTS:
    sn_reps = get_replicates(tpm_df, construct, "SN")
    lys_reps = get_replicates(tpm_df, construct, "Lysate")

    averaged_data[f"SN_{construct}_avg"] = sn_reps.mean(axis=1)
    averaged_data[f"Lysate_{construct}_avg"] = lys_reps.mean(axis=1)


gtf_data = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_is_mt = gtf_data["gene_id_is_mt"]
averaged_data["is_mt"] = [gene_id_is_mt.get(gene_id) for gene_id in averaged_data.index]

# Use GridSpec: 5 equal subplot columns + 1 narrow colorbar column
fig = plt.figure(figsize=(14, 2.5))
gs = GridSpec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.15], wspace=0.3)
axes = [fig.add_subplot(gs[0, i]) for i in range(5)]
cbar_ax = fig.add_subplot(gs[0, 5])

last_hb = None
for i, construct in enumerate(CONSTRUCTS):
    sn_col = f"SN_{construct}_avg"
    lys_col = f"Lysate_{construct}_avg"
    title = NAME_MAPPING[construct]
    last_hb = plot_hexbin(
        axes[i],
        averaged_data[lys_col],
        averaged_data[sn_col],
        title,
        is_mt=averaged_data["is_mt"],
    )

# Colorbar in its own dedicated axes
fig.colorbar(last_hb, cax=cbar_ax, label="Count (log)")

# Legend with MT marker on first subplot only
handles = [
    plt.Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor="red",
        markeredgecolor="darkred",
        markersize=8,
        label="MT genes",
    )
]
axes[0].legend(handles=handles, loc="lower right", fontsize=PLOT_FONT_SIZE)

save(fig, "suppl11b_averaged_hexbin")
plt.close(fig)
print("Done.")
