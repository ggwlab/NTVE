"""
Supplementary Figure 21b,c upstream rankings.

Recreates the ranking table from:
- Suppl21/Rank_Genes_Per_Timepoint.ipynb

Input:
- Suppl21/cardio_deseq2_cubicsplines_output/fitted_trajectories_all_genes.csv

Outputs:
- refactoring_roadmap/suppl21bc_output/all_gene_rankings_log2_centered.csv

Notebook helper output remains available in code, but is disabled by default:
- refactoring_roadmap/suppl21bc_plots/top_ranked_genes_heatmap_all.{png,svg}
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def find_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here.parent, here.parent.parent):
        if (candidate / "resources").exists():
            return candidate
    return here.parent.parent

ROOT = find_root()
TRAJ_CSV = ROOT / "Figure7" / "cardio_deseq2_cubicsplines_output" / "fitted_trajectories_all_genes.csv"
OUT_DATA = Path(__file__).parent / "suppl21bc_output"
OUT_PLOTS = Path(__file__).parent / "suppl21bc_plots"
OUT_DATA.mkdir(exist_ok=True)
OUT_PLOTS.mkdir(exist_ok=True)

PSEUDOCOUNT = 0.1


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "svg"):
        out = OUT_PLOTS / f"{stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


def create_rankings_for_timepoint_all_genes(timepoint: int, obs_log_centered: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    values = obs_log_centered[timepoint].copy()
    upregulated = values[values > 0].sort_values(ascending=False)
    downregulated = values[values < 0].sort_values(ascending=True)

    up_df = pd.DataFrame(
        {
            "Gene": upregulated.index,
            "Log2-Centered": upregulated.values,
            "Rank": range(1, len(upregulated) + 1),
            "Direction": "Upregulated",
        }
    ).reset_index(drop=True)
    down_df = pd.DataFrame(
        {
            "Gene": downregulated.index,
            "Log2-Centered": downregulated.values,
            "Rank": range(1, len(downregulated) + 1),
            "Direction": "Downregulated",
        }
    ).reset_index(drop=True)
    return up_df, down_df


print("Loading fitted trajectories...")
if not TRAJ_CSV.exists():
    raise FileNotFoundError(
        "Missing fitted cardio trajectories. Run `fig7c_pipeline.py` first to generate "
        f"{ROOT / 'Figure7' / 'cardio_deseq2_cubicsplines_output' / 'fitted_trajectories_all_genes.csv'}."
    )
trajectories = pd.read_csv(TRAJ_CSV)
print(f"  Total genes: {trajectories['Gene'].nunique()}")
print(f"  Timepoints: {sorted(trajectories['Time'].unique())}")

heatmap_data_obs = trajectories.pivot(index="Gene", columns="Time", values="Observed")
obs_log = np.log2(heatmap_data_obs + PSEUDOCOUNT)
obs_log_centered = obs_log.sub(obs_log.mean(axis=1), axis=0)

timepoints = sorted(obs_log_centered.columns)
all_rankings: dict[int, dict[str, pd.DataFrame]] = {}
for tp in timepoints:
    up, down = create_rankings_for_timepoint_all_genes(tp, obs_log_centered)
    all_rankings[tp] = {"upregulated": up, "downregulated": down}
    print(f"Day {tp}: {len(up)} upregulated, {len(down)} downregulated genes")

all_data = []
for tp in timepoints:
    up_df = all_rankings[tp]["upregulated"].copy()
    down_df = all_rankings[tp]["downregulated"].copy()
    up_df["Timepoint"] = tp
    down_df["Timepoint"] = tp
    all_data.append(up_df)
    all_data.append(down_df)

combined_rankings = pd.concat(all_data, ignore_index=True)
combined_rankings = combined_rankings[["Timepoint", "Direction", "Rank", "Gene", "Log2-Centered"]]
rankings_csv = OUT_DATA / "all_gene_rankings_log2_centered.csv"
combined_rankings.to_csv(rankings_csv, index=False)
print(f"Saved: {rankings_csv}")

# Notebook helper output intentionally disabled by default:
#
# top_genes_per_tp = {}
# for tp in timepoints:
#     up = all_rankings[tp]["upregulated"].head(20)["Gene"].tolist()
#     down = all_rankings[tp]["downregulated"].head(20)["Gene"].tolist()
#     top_genes_per_tp[tp] = up + down
#
# all_top_genes = set()
# for genes in top_genes_per_tp.values():
#     all_top_genes.update(genes)
#
# heatmap_data = obs_log_centered.loc[obs_log_centered.index.isin(list(all_top_genes))]
# sort_order = heatmap_data.idxmax(axis=1).sort_values().index
# heatmap_data = heatmap_data.loc[sort_order]
#
# fig, ax = plt.subplots(figsize=(12, max(8, len(heatmap_data) * 0.15)))
# im = ax.imshow(
#     heatmap_data.values,
#     cmap="RdBu_r",
#     aspect="auto",
#     vmin=np.percentile(heatmap_data.values, 2),
#     vmax=np.percentile(heatmap_data.values, 98),
#     interpolation="nearest",
# )
# ax.set_xlabel("Timepoint (days)", fontweight="bold")
# ax.set_ylabel("Gene", fontweight="bold")
# ax.set_title("Top Ranked Genes per Timepoint - Cardio Development\n(Log2-Centered, Sorted by Peak Time)", fontweight="bold")
# ax.set_xticks(range(len(heatmap_data.columns)))
# ax.set_xticklabels(heatmap_data.columns.astype(int), rotation=0)
# ax.set_yticks(range(len(heatmap_data)))
# ax.set_yticklabels(heatmap_data.index, fontsize=min(8, 100 / len(heatmap_data)))
# plt.colorbar(im, ax=ax, label="log₂(Expression) - Mean")
# fig.tight_layout()
# save(fig, "top_ranked_genes_heatmap_all")
# plt.close(fig)

print("Done.")
