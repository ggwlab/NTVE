"""
Figure 2a — Benchmarking comparison total RNA amount in supernatant
Standalone reproduction of the mGL-based RNA amount panel from
Suppl11/benchmarking_simplified_final.ipynb.
Outputs to: refactoring_roadmap/2a_plots/
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir

DATA_DIR = get_generated_figure2_loaded_dir()
OUT_DIR = Path(__file__).parent / "2a_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "2a_csv"
CSV_DIR.mkdir(exist_ok=True)

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


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


print("Loading benchmarking input estimates...")
tpm_df = pd.read_csv(DATA_DIR / "elo_tpm_matrix_protein_coding.csv", index_col=0)
gene_names = tpm_df["GeneName"]
numreads_df = pd.read_csv(DATA_DIR / "elo_numreads_matrix_protein_coding.csv", index_col=0)

mgl_reads = numreads_df[gene_names == "mGreenLantern"]
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
log_stats_df = (
    mgl_plot_df.assign(log10_Input_ng=np.log10(mgl_plot_df["Input_ng"]))
    .groupby("Condition")["log10_Input_ng"]
    .agg(["mean", "std"])
    .reset_index()
)
stats_df = log_stats_df.rename(columns={"mean": "log10_mean", "std": "log10_std"})
stats_df["mean"] = 10 ** stats_df["log10_mean"]
stats_df["lower"] = 10 ** (stats_df["log10_mean"] - stats_df["log10_std"])
stats_df["upper"] = 10 ** (stats_df["log10_mean"] + stats_df["log10_std"])
stats_df["yerr_lower"] = stats_df["mean"] - stats_df["lower"]
stats_df["yerr_upper"] = stats_df["upper"] - stats_df["mean"]
stats_df["Condition"] = pd.Categorical(
    stats_df["Condition"],
    categories=[NAME_MAPPING[c] for c in CONSTRUCTS],
    ordered=True,
)
stats_df = stats_df.sort_values("Condition")
stats_df["n_replicates"] = stats_df["Condition"].map(mgl_plot_df.groupby("Condition").size())
mgl_plot_df.to_csv(CSV_DIR / "2a_total_rna_amount_points.csv", index=False)
stats_df.to_csv(CSV_DIR / "2a_total_rna_amount_summary.csv", index=False)

fig, ax = plt.subplots(figsize=(PLOT_WIDTH * 0.75, PLOT_HEIGHT))
x_pos = np.arange(len(stats_df))
bar_width = 0.55
sn_color = "#3498db"

ax.bar(
    x_pos,
    stats_df["mean"],
    bar_width,
    yerr=np.vstack([stats_df["yerr_lower"], stats_df["yerr_upper"]]),
    capsize=5,
    color=sn_color,
    alpha=0.7,
    error_kw={"linewidth": 1.5},
)

for i, condition in enumerate(stats_df["Condition"]):
    pts = mgl_plot_df[mgl_plot_df["Condition"] == condition]["Input_ng"]
    jitter = np.random.normal(0, 0.05, size=len(pts))
    ax.scatter(np.ones(len(pts)) * i + jitter, pts, color="black", s=15, alpha=0.6, zorder=10)

ax.set_yscale("log")
ax.axhline(1, ls="--", color="gray", alpha=0.5, linewidth=1.5, label="1 ng reference")
ax.set_title(
    "Benchmarking Experiment - Estimated RNA Input (SN only)\n(Bars: Geometric Mean ± 1 log10-SD, Points: Individual Replicates)",
    fontweight="bold",
    fontsize=PLOT_FONT_SIZE + 1,
)
ax.set_ylabel("Estimated RNA Input (ng)", fontweight="bold")
ax.set_xlabel("Condition", fontweight="bold")
ax.set_xticks(x_pos)
ax.set_xticklabels(stats_df["Condition"], rotation=45, ha="right")
ax.legend(fontsize=PLOT_FONT_SIZE)
ax.grid(True, alpha=0.3, axis="y")
fig.tight_layout()

save(fig, "2a_total_rna_amount_barplot")
plt.close(fig)
print("Done.")
