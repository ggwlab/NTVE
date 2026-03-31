"""
Supplementary Figure 10c-e — STAR vs Minimap benchmarking comparison
Standalone reproduction of the workbook-linked outputs from
Suppl10/compare_STAR_vs_Minimap.ipynb.

Workbook mapping verified from SupplementaryTable6_Numeric data.xlsx:
  - Suppl 10c   -> Suppl10/plots/Venn_detected_genes.svg
  - Suppl 10d,e -> Suppl10/plots/3D_replicate_comparisons.svg

This script also reproduces the notebook CSV side export:
  - csv/STAR_detected_genes.csv
  - csv/Minimap_detected_genes.csv
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn2
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


ROOT = Path(__file__).parent.parent

from ntvetools import load_gtf_df


def first_existing_dir(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


DATA_DIR = first_existing_dir(ROOT / "resources" / "sup5_comparison_data", ROOT / "Suppl10" / "comparison_data")
OUT_DIR = Path(__file__).parent / "suppl10_plots"
CSV_DIR = Path(__file__).parent / "suppl10_csv"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR.mkdir(exist_ok=True)

plt.rcParams["svg.fonttype"] = "none"

STAR_LYS_COLS = ["24L006460", "24L006461", "24L006462"]
STAR_SN_COLS = ["24L006479", "24L006480", "24L006481"]
MINIMAP_LYS_COLS = ["24L006460", "24L006461", "24L006462"]
MINIMAP_SN_COLS = ["24L006479", "24L006480", "24L006481"]


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=300)
        print(f"Saved: {out}")


def indices_gt_zero(series: pd.Series) -> set[str]:
    return set(series[series.gt(0)].index)


def create_detection_dataframe(
    lys_detected: set[str],
    sn_detected: set[str],
    gene_id_to_name: dict[str, str],
) -> pd.DataFrame:
    lys_only = sorted(lys_detected - sn_detected)
    shared = sorted(lys_detected & sn_detected)
    ntve_only = sorted(sn_detected - lys_detected)

    rows = []
    for gene in lys_only:
        rows.append({"ENSG_id": gene, "Gene_Name": gene_id_to_name.get(gene, gene), "Category": "lysate_only"})
    for gene in shared:
        rows.append({"ENSG_id": gene, "Gene_Name": gene_id_to_name.get(gene, gene), "Category": "shared"})
    for gene in ntve_only:
        rows.append({"ENSG_id": gene, "Gene_Name": gene_id_to_name.get(gene, gene), "Category": "ntve_only"})
    return pd.DataFrame(rows)


print("Loading Suppl10 comparison inputs...")
star_tpm = pd.read_csv(DATA_DIR / "STAR_protein_coding_tpm_normalized.csv", index_col=0)[STAR_LYS_COLS + STAR_SN_COLS]
minimap_tpm = pd.read_csv(DATA_DIR / "Minimap_protein_coding_tpm_normalized.csv", index_col=0)
common_genes = [gene for gene in minimap_tpm.index if gene in star_tpm.index]
star_tpm_common = star_tpm.loc[common_genes]
minimap_tpm_common = minimap_tpm.loc[common_genes]

print("Loading gene-name annotation for CSV exports...")
gtf_data = load_gtf_df(str(ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_name = gtf_data["gene_id_to_name"]


print("Reproducing Suppl 10c: detected-gene Venn diagram...")
star_lys_detected = set(star_tpm_common.index[(star_tpm_common[STAR_LYS_COLS] > 0).all(axis=1)])
star_sn_detected = set(star_tpm_common.index[(star_tpm_common[STAR_SN_COLS] > 0).all(axis=1)])
minimap_lys_detected = set(minimap_tpm_common.index[(minimap_tpm_common[MINIMAP_LYS_COLS] > 0).all(axis=1)])
minimap_sn_detected = set(minimap_tpm_common.index[(minimap_tpm_common[MINIMAP_SN_COLS] > 0).all(axis=1)])

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

venn2(
    [minimap_lys_detected, minimap_sn_detected],
    set_labels=("Lysate", "NTVE"),
    ax=axes[0],
)
axes[0].set_title(
    "Minimap: Genes Detected in All Replicates\n(TPM > 0 in all 3 replicates)",
    fontsize=14,
    fontweight="bold",
)

venn2(
    [star_lys_detected, star_sn_detected],
    set_labels=("Lysate", "NTVE"),
    ax=axes[1],
)
axes[1].set_title(
    "STAR: Genes Detected in All Replicates\n(TPM > 0 in all 3 replicates)",
    fontsize=14,
    fontweight="bold",
)

plt.tight_layout()
save(fig, "Venn_detected_genes")
plt.close(fig)


print("Writing notebook-faithful detected-gene CSV tables...")
star_df = create_detection_dataframe(star_lys_detected, star_sn_detected, gene_id_to_name)
star_csv = CSV_DIR / "STAR_detected_genes.csv"
star_df.to_csv(star_csv, index=False)
print(f"Saved: {star_csv}")

minimap_df = create_detection_dataframe(minimap_lys_detected, minimap_sn_detected, gene_id_to_name)
minimap_csv = CSV_DIR / "Minimap_detected_genes.csv"
minimap_df.to_csv(minimap_csv, index=False)
print(f"Saved: {minimap_csv}")


print("Reproducing Suppl 10d,e: 3D replicate comparisons...")
star_lys_data = star_tpm_common[STAR_LYS_COLS]
star_sn_data = star_tpm_common[STAR_SN_COLS]
minimap_lys_data = minimap_tpm_common[MINIMAP_LYS_COLS]
minimap_sn_data = minimap_tpm_common[MINIMAP_SN_COLS]

log10_star_lys = np.log10(star_lys_data.values + 1)
log10_star_sn = np.log10(star_sn_data.values + 1)
log10_minimap_lys = np.log10(minimap_lys_data.values + 1)
log10_minimap_sn = np.log10(minimap_sn_data.values + 1)

all_data = np.concatenate(
    [
        log10_star_lys.ravel(),
        log10_star_sn.ravel(),
        log10_minimap_lys.ravel(),
        log10_minimap_sn.ravel(),
    ]
)
global_min = np.floor(all_data.min())
global_max = np.ceil(all_data.max())
ticks = np.arange(global_min, global_max + 1, 1)

fig = plt.figure(figsize=(8, 7))

panels = [
    (221, log10_star_lys, "STAR Lysate Replicates", ("Lysate Rep. 1", "Lysate Rep. 2", "Lysate Rep. 3")),
    (222, log10_star_sn, "STAR NTVE Replicates", ("NTVE Rep. 1", "NTVE Rep. 2", "NTVE Rep. 3")),
    (223, log10_minimap_lys, "Minimap Lysate Replicates", ("Lysate Rep. 1", "Lysate Rep. 2", "Lysate Rep. 3")),
    (224, log10_minimap_sn, "Minimap NTVE Replicates", ("NTVE Rep. 1", "NTVE Rep. 2", "NTVE Rep. 3")),
]

for subplot_code, values, title, labels in panels:
    ax = fig.add_subplot(subplot_code, projection="3d")
    color_values = np.mean(values, axis=1)
    scatter = ax.scatter(
        values[:, 0],
        values[:, 1],
        values[:, 2],
        c=color_values,
        cmap="inferno",
        s=1,
        alpha=0.6,
        rasterized=True,
    )
    ax.set_xlabel(labels[0], fontsize=8)
    ax.set_ylabel(labels[1], fontsize=8)
    ax.set_zlabel(labels[2], fontsize=8)
    ax.set_title(title, fontsize=8, fontweight="bold")
    fig.colorbar(scatter, ax=ax, label="Mean log10(TPM+1)", shrink=0.5, pad=0.1).ax.tick_params(labelsize=8)
    ax.set_xlim(global_min, global_max)
    ax.set_ylim(global_min, global_max)
    ax.set_zlim(global_min, global_max)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_zticks(ticks)
    ax.tick_params(axis="both", which="major", labelsize=8)

fig.suptitle("3D Replicate Comparisons: Each Axis = One Replicate", fontsize=8, fontweight="bold")
plt.tight_layout(w_pad=2.0, h_pad=2.0)
save(fig, "3D_replicate_comparisons")
plt.close(fig)

print("Done.")
