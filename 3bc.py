"""
Figure 3b,c — ROC and Precision-Recall curves
Notebook-faithful analysis with manuscript-facing combined export.
Outputs to: refactoring_roadmap/3bc_plots/3bc_roc_pr.svg and .png
"""

from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.utils import resample

ROOT = Path(__file__).parent.parent


def first_existing_dir(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


FIG3 = first_existing_dir(ROOT / "resources" / "fig3", ROOT / "Figure3")
OUT_DIR = Path(__file__).parent / "3bc_plots"
OUT_DIR.mkdir(exist_ok=True)


def load_data(file_name: str) -> pd.DataFrame:
    return pd.read_csv(FIG3 / file_name).rename(columns={"Unnamed: 0": "GeneID", "Geneid": "GeneID"}).dropna()


def create_label_set(df: pd.DataFrame, threshold: float = 0.05, top_n: int = 20) -> set[str]:
    top_genes = df.sort_values("padj").head(top_n)
    gene_set = set(top_genes["GeneID"])
    threshold_set = set(df.query(f"padj < {threshold}")["GeneID"])
    if len(threshold_set) > top_n:
        return threshold_set
    return gene_set


def downsample_majority(positive_set: set[str], negative_set: set[str]) -> tuple[set[str], set[str]]:
    positive_list = list(positive_set)
    negative_list = list(negative_set)
    if len(negative_list) > len(positive_list):
        majority_list = negative_list
        minority_list = positive_list
    else:
        majority_list = positive_list
        minority_list = negative_list

    majority_downsampled = resample(
        majority_list,
        replace=False,
        n_samples=len(minority_list),
        random_state=42,
    )

    if len(negative_list) > len(positive_list):
        return set(positive_list), set(majority_downsampled)
    return set(majority_downsampled), set(minority_list)


def calculate_curves(positive_set: set[str], negative_set: set[str], predictor_df: pd.DataFrame, downsample: bool = False):
    if downsample:
        positive_set, negative_set = downsample_majority(positive_set, negative_set)

    all_genes = list(positive_set.union(negative_set))
    labels = [1 if gene in positive_set else 0 for gene in all_genes]

    merged_df = predictor_df[predictor_df["GeneID"].isin(all_genes)].copy()
    merged_df["label"] = merged_df["GeneID"].map(dict(zip(all_genes, labels)))
    merged_df["prediction"] = abs(merged_df["log2FoldChange"])
    merged_df = merged_df.dropna(subset=["prediction", "label"])

    if len(merged_df) == 0 or len(set(merged_df["label"])) < 2:
        return None, None, None, None, None, None, None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fpr, tpr, _ = roc_curve(merged_df["label"], merged_df["prediction"])
        roc_auc = auc(fpr, tpr)
        precision, recall, _ = precision_recall_curve(merged_df["label"], merged_df["prediction"])
        pr_auc = average_precision_score(merged_df["label"], merged_df["prediction"])

    return fpr, tpr, roc_auc, precision, recall, pr_auc, merged_df["label"].values


def plot_combined(curve_data: dict, save_stem: str) -> None:
    condition = list(curve_data.keys())[0]
    data = curve_data[condition]
    if data[0] is None:
        raise RuntimeError("No valid ROC/PR data to plot")

    fpr, tpr, roc_auc, precision, recall, pr_auc, labels = data
    baseline = np.mean(labels)

    x_axes_size_mm = 50
    y_axes_size_mm = 25
    margin_mm = 10
    panel_gap_mm = 8
    figure_width_mm = x_axes_size_mm * 2 + 2 * margin_mm + panel_gap_mm
    figure_height_mm = y_axes_size_mm + 2 * margin_mm
    fig = plt.figure(figsize=(figure_width_mm / 25.4, figure_height_mm / 25.4), dpi=300)

    plt.rc("font", size=8)
    plt.rc("axes", titlesize=8)
    plt.rc("axes", labelsize=8)
    plt.rc("xtick", labelsize=8)
    plt.rc("ytick", labelsize=8)
    plt.rc("legend", fontsize=6)
    plt.rcParams["svg.fonttype"] = "none"

    left = margin_mm / figure_width_mm
    bottom = margin_mm / figure_height_mm
    width = x_axes_size_mm / figure_width_mm
    height = y_axes_size_mm / figure_height_mm
    gap = panel_gap_mm / figure_width_mm

    ax1 = fig.add_axes([left, bottom, width, height])
    ax2 = fig.add_axes([left + width + gap, bottom, width, height])

    ax1.plot(fpr, tpr, lw=1.5, color="#1f77b4", label=f"{condition} (AUC = {roc_auc:.2f})")
    ax1.plot([0, 1], [0, 1], color="gray", lw=1, linestyle="--", label="Random Classifier")
    ax1.set_xlabel("False Positive Rate")
    ax1.set_ylabel("True Positive Rate")
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.legend(loc="lower right")

    ax2.plot(recall, precision, lw=1.5, color="#ff7f0e", label=f"{condition} (AP = {pr_auc:.2f})")
    ax2.axhline(y=baseline, color="#ff7f0e", linestyle="--", alpha=0.7, lw=1, label=f"Random (y={baseline:.2f})")
    ax2.set_xlabel("Recall")
    ax2.set_ylabel("Precision")
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.legend(loc="upper right")

    for ext in ("svg", "png"):
        out = OUT_DIR / f"{save_stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


conditions = ["IFNg"]
curve_data = {}
curve_data_downsampled = {}

for condition in conditions:
    lysate_gt_df = load_data(f"greaterAbs_Dox_vs_{condition}_lysate_base_mean_gt_0.csv")
    lysate_lt_df = load_data(f"lessAbs_Dox_vs_{condition}_lysate_base_mean_gt_0.csv")
    positive_set = create_label_set(lysate_gt_df, top_n=39)
    negative_set = create_label_set(lysate_lt_df, top_n=39)
    sn_df = load_data(f"greaterAbs_Dox_vs_{condition}_SN_base_mean_gt_0.csv")

    curve_data[condition] = calculate_curves(positive_set, negative_set, sn_df, downsample=False)
    curve_data_downsampled[condition] = calculate_curves(positive_set, negative_set, sn_df, downsample=True)

    print(f"\nCondition: {condition}")
    print(f"Original - Positive: {len(positive_set)}, Negative: {len(negative_set)}")
    pos_down, neg_down = downsample_majority(positive_set, negative_set)
    print(f"Downsampled - Positive: {len(pos_down)}, Negative: {len(neg_down)}")

print("\nAUC Values:")
for condition in conditions:
    original = curve_data[condition]
    down = curve_data_downsampled[condition]
    print(f"\n{condition}:")
    print(f"  Original    - ROC AUC = {original[2]:.4f}, PR AUC = {original[5]:.4f}")
    print(f"  Downsampled - ROC AUC = {down[2]:.4f}, PR AUC = {down[5]:.4f}")

plot_combined(curve_data, "3bc_roc_pr")
print("Done.")
