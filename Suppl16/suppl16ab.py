"""
Supplementary Figure 16a, 16b — FLAG purification correlation heatmaps.

Source notebook: Suppl16/Purification_Comparison-simplified.ipynb

Manuscript-faithful default outputs:
  suppl16a_flag_purification_correlation  — human Pearson 6×6 SN panel
  suppl16b_lysate_vs_ntve_correlation     — human Pearson 3×6 Lysate vs SN panel

Additional notebook helper heatmaps remain available in functions below, but are
not emitted by default.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
RESOURCES = REPO_ROOT / "resources"
OUT_DIR = Path(__file__).resolve().parent / "suppl16ab_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).resolve().parent / "suppl16ab_csv"
CSV_DIR.mkdir(exist_ok=True)

plt.rcParams["svg.fonttype"] = "none"

# ── Load featureCounts files ───────────────────────────────────────────────────
fc1 = pd.read_csv(
    RESOURCES / "harmonized_harmonized_gene_counts_rv_stranded.txt",
    sep="\t", skiprows=1,
).set_index("Geneid")

fc2 = pd.read_csv(
    RESOURCES / "comparison_harmonized_gene_counts.txt.gz",
    sep="\t", skiprows=1,
).set_index("Geneid")

sample_cols1 = sorted(fc1.columns[5:])
sample_cols2 = sorted(fc2.columns[5:])

raw1 = fc1.rename(columns={c: c.split("/")[-1].split("_")[0] for c in sample_cols1})
raw2 = fc2.rename(columns={
    c: c.split("/")[-1].replace("_untrimmed_harmonizedAligned.sortedByCoord.out.bam", "")
    for c in sample_cols2
})

# Merge: keep metadata cols from fc1, sample cols from both
raw_counts = pd.concat(
    [raw1.iloc[:, :5], raw1.iloc[:, 5:], raw2.iloc[:, 5:]],
    axis=1,
)

# ── Filter to gene-level (human ENSG* or mouse MUS*) ──────────────────────────
gene_df = raw_counts.loc[
    raw_counts.index.str.contains("MUS") | raw_counts.index.str.contains("ENSG")
]

# All sample IDs used across both experiments
all_sample_ids = [
    "24L006460", "24L006461", "24L006462",   # pure Lysate HEK
    "24L006479", "24L006480", "24L006481",   # pure SN HEK
    "24L010907", "24L010908", "24L010909",   # SN N2a
    "24L010910", "24L010911", "24L010912",   # SN HEK
    "Pulldown_S1_1", "Pulldown_S1_2", "Pulldown_S1_3",  # Lysate HEK
    "Pulldown_S2_1", "Pulldown_S2_2", "Pulldown_S2_3",  # Lysate N2a
    "Pulldown_S3_1", "Pulldown_S3_2", "Pulldown_S3_3",  # Lysate co-culture
    "Pulldown_S4_1", "Pulldown_S4_2", "Pulldown_S4_3",  # Lysate DMSO
]

# ── Compute RPM (library-size normalise across all samples) ────────────────────
rpm_df = gene_df[all_sample_ids] / (gene_df[all_sample_ids].sum(axis=0) * 1e-6)

rpm_homo = rpm_df.loc[rpm_df.index.str.contains("ENSG")]
rpm_mus  = rpm_df.loc[rpm_df.index.str.contains("MUS")]

# ── Renormalise per-species so each sample sums to 1 M ────────────────────────
def renormalize_rpm(df: pd.DataFrame) -> pd.DataFrame:
    return df * (1_000_000 / df.sum(axis=0))

rpm_homo_norm = renormalize_rpm(rpm_homo)
rpm_mus_norm  = renormalize_rpm(rpm_mus)


# ── Helper: save figure as svg + png ──────────────────────────────────────────
def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {stem}.svg / .png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 — sn_lysate  (Pearson, HEK + N2a only, no pure refs)
# ══════════════════════════════════════════════════════════════════════════════
def plot_sn_lysate(df: pd.DataFrame, species: str) -> plt.Figure:
    organized = {
        "SN": {
            "HEK": ["24L010910", "24L010911", "24L010912"],
            "N2a": ["24L010907", "24L010908", "24L010909"],
        },
        "Lysate": {
            "HEK":        ["Pulldown_S1_1", "Pulldown_S1_2", "Pulldown_S1_3"],
            "N2a":        ["Pulldown_S2_1", "Pulldown_S2_2", "Pulldown_S2_3"],
            "co-culture": ["Pulldown_S3_1", "Pulldown_S3_2", "Pulldown_S3_3"],
            "DMSO":       ["Pulldown_S4_1", "Pulldown_S4_2", "Pulldown_S4_3"],
        },
    }
    condition_order = ["HEK", "N2a", "co-culture", "DMSO"]
    return _heatmap_conditions(
        df, organized, condition_order, species,
        method="pearson", figsize=(16, 14),
        title_suffix="(SN vs Lysate Samples by Condition)",
    )


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 — all_samples  (Spearman, all samples incl. pure HEK refs)
# ══════════════════════════════════════════════════════════════════════════════
def plot_all_samples(df: pd.DataFrame, species: str) -> plt.Figure:
    organized = {
        "SN": {
            "pure_HEK": ["24L006479", "24L006480", "24L006481"],
            "HEK":      ["24L010910", "24L010911", "24L010912"],
            "N2a":      ["24L010907", "24L010908", "24L010909"],
        },
        "Lysate": {
            "pure_HEK":   ["24L006460", "24L006461", "24L006462"],
            "HEK":        ["Pulldown_S1_1", "Pulldown_S1_2", "Pulldown_S1_3"],
            "N2a":        ["Pulldown_S2_1", "Pulldown_S2_2", "Pulldown_S2_3"],
            "co-culture": ["Pulldown_S3_1", "Pulldown_S3_2", "Pulldown_S3_3"],
            "DMSO":       ["Pulldown_S4_1", "Pulldown_S4_2", "Pulldown_S4_3"],
        },
    }
    condition_order = ["pure_HEK", "HEK", "N2a", "co-culture", "DMSO"]
    return _heatmap_conditions(
        df, organized, condition_order, species,
        method="spearman", figsize=(18, 16),
        title_suffix="(All Samples incl. Pure HEK refs)",
    )


def _heatmap_conditions(
    df, organized, condition_order, species, method, figsize, title_suffix
) -> plt.Figure:
    ordered_cols, labels, boundaries = [], [], []
    pos = 0

    for sample_type in ["SN", "Lysate"]:
        for cond in condition_order:
            if cond not in organized[sample_type]:
                continue
            existing = [s for s in organized[sample_type][cond] if s in df.columns]
            if not existing:
                continue
            for i, s in enumerate(existing):
                ordered_cols.append(s)
                labels.append(f"{sample_type}_{cond}\nRep{i+1}")
            pos += len(existing)
        boundaries.append(pos)

    if not ordered_cols:
        raise ValueError(f"No matching columns found for {species}")

    corr = df[ordered_cols].corr(method=method)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        corr, annot=True, cmap="inferno", vmin=0.5, vmax=1.0, fmt=".2f",
        xticklabels=labels, yticklabels=labels, square=True, ax=ax,
    )

    # white separators between condition groups
    for b in boundaries[:-1]:
        ax.axhline(y=b, color="white", linewidth=3)
        ax.axvline(x=b, color="white", linewidth=3)

    # yellow line between SN block and Lysate block
    sn_end = sum(
        len([s for s in organized["SN"][c] if s in df.columns])
        for c in condition_order if c in organized["SN"]
    )
    if sn_end:
        ax.axhline(y=sn_end, color="yellow", linewidth=4)
        ax.axvline(x=sn_end, color="yellow", linewidth=4)

    ax.set_title(
        f"{species} Gene Expression Correlation Heatmap\n{title_suffix}",
        fontsize=16, pad=20,
    )
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3 — 6×6 SN  (Pearson, pure_HEK SN vs HEK SN)
# ══════════════════════════════════════════════════════════════════════════════
def plot_6x6_sn(df: pd.DataFrame, species: str) -> plt.Figure:
    cols = [
        "24L006479", "24L006480", "24L006481",   # pure_HEK SN
        "24L010910", "24L010911", "24L010912",   # HEK SN
    ]
    labels = [
        "pure_SN_HEK\nRep1", "pure_SN_HEK\nRep2", "pure_SN_HEK\nRep3",
        "SN_HEK\nRep1",      "SN_HEK\nRep2",      "SN_HEK\nRep3",
    ]
    corr = df[cols].corr(method="pearson")

    fig, ax = plt.subplots(figsize=(4, 4))
    sns.heatmap(
        corr, annot=True, cmap="inferno", vmin=0.0, vmax=1.0, fmt=".2f",
        xticklabels=labels, yticklabels=labels, square=True, ax=ax,
        annot_kws={"size": 8}, cbar_kws={"shrink": 0.8},
    )
    ax.axhline(y=3, color="white", linewidth=2)
    ax.axvline(x=3, color="white", linewidth=2)
    ax.set_title(
        f"{species} – SN Correlation (6×6)\npure_HEK vs HEK", fontsize=8, pad=8
    )
    ax.tick_params(labelsize=8)
    ax.figure.axes[1].tick_params(labelsize=8)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 — 3×6 Lysate vs SN  (Pearson, pure_HEK Lysate vs SN groups)
# ══════════════════════════════════════════════════════════════════════════════
def plot_3x6_lys_sn(df: pd.DataFrame, species: str) -> plt.Figure:
    lysate_3 = ["24L006460", "24L006461", "24L006462"]
    sn_6 = [
        "24L006479", "24L006480", "24L006481",
        "24L010910", "24L010911", "24L010912",
    ]
    row_labels = ["pure_Lys_HEK\nRep1", "pure_Lys_HEK\nRep2", "pure_Lys_HEK\nRep3"]
    col_labels = [
        "pure_SN_HEK\nRep1", "pure_SN_HEK\nRep2", "pure_SN_HEK\nRep3",
        "SN_HEK\nRep1",      "SN_HEK\nRep2",      "SN_HEK\nRep3",
    ]
    cross_corr = df[lysate_3 + sn_6].corr(method="pearson").loc[lysate_3, sn_6]

    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(
        cross_corr, annot=True, cmap="inferno", vmin=0.0, vmax=1.0, fmt=".2f",
        xticklabels=col_labels, yticklabels=row_labels, ax=ax,
        annot_kws={"size": 8}, cbar_kws={"shrink": 0.8},
    )
    ax.axvline(x=3, color="white", linewidth=2)
    ax.set_title(
        f"{species} – Lysate vs SN Correlation (3×6)\npure_HEK Lysate vs SN",
        fontsize=8, pad=8,
    )
    ax.tick_params(labelsize=8)
    ax.figure.axes[1].tick_params(labelsize=8)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    fig.tight_layout()
    return fig


# ── Main ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n── Human manuscript panels ──")

    # 16a: 6×6 SN correlation matrix
    cols_6x6 = [
        "24L006479", "24L006480", "24L006481",
        "24L010910", "24L010911", "24L010912",
    ]
    labels_6x6 = [
        "pure_SN_HEK_Rep1", "pure_SN_HEK_Rep2", "pure_SN_HEK_Rep3",
        "SN_HEK_Rep1", "SN_HEK_Rep2", "SN_HEK_Rep3",
    ]
    corr_6x6 = rpm_homo_norm[cols_6x6].corr(method="pearson")
    corr_6x6.index = labels_6x6
    corr_6x6.columns = labels_6x6
    corr_6x6.index.name = "sample"
    corr_6x6.to_csv(CSV_DIR / "suppl16a_pearson_6x6_sn.csv")
    print(f"  saved suppl16a_pearson_6x6_sn.csv")

    # 16b: 3×6 Lysate vs SN cross-correlation matrix
    lysate_3 = ["24L006460", "24L006461", "24L006462"]
    sn_6 = ["24L006479", "24L006480", "24L006481", "24L010910", "24L010911", "24L010912"]
    row_labels_3x6 = ["pure_Lys_HEK_Rep1", "pure_Lys_HEK_Rep2", "pure_Lys_HEK_Rep3"]
    col_labels_3x6 = [
        "pure_SN_HEK_Rep1", "pure_SN_HEK_Rep2", "pure_SN_HEK_Rep3",
        "SN_HEK_Rep1", "SN_HEK_Rep2", "SN_HEK_Rep3",
    ]
    cross_corr = rpm_homo_norm[lysate_3 + sn_6].corr(method="pearson").loc[lysate_3, sn_6]
    cross_corr.index = row_labels_3x6
    cross_corr.columns = col_labels_3x6
    cross_corr.to_csv(CSV_DIR / "suppl16b_pearson_3x6_lysate_vs_sn.csv")
    print(f"  saved suppl16b_pearson_3x6_lysate_vs_sn.csv")

    save(
        plot_6x6_sn(rpm_homo_norm, "Human"),
        "suppl16a_flag_purification_correlation",
    )
    save(
        plot_3x6_lys_sn(rpm_homo_norm, "Human"),
        "suppl16b_lysate_vs_ntve_correlation",
    )

    # Notebook helper outputs are intentionally disabled by default so the
    # refactoring output directory reflects the manuscript panel set.
    #
    # save(plot_sn_lysate(rpm_homo_norm, "Human"), "correlation_heatmap_human_sn_lysate")
    # save(plot_all_samples(rpm_homo_norm, "Human"), "correlation_heatmap_human_all_samples")
    # save(plot_6x6_sn(rpm_homo_norm, "Human"), "correlation_heatmap_human_6x6_sn")
    # save(plot_3x6_lys_sn(rpm_homo_norm, "Human"), "correlation_heatmap_human_3x6_lys_sn")
    # save(plot_sn_lysate(rpm_mus_norm, "Mouse"), "correlation_heatmap_mouse_sn_lysate")
    # save(plot_all_samples(rpm_mus_norm, "Mouse"), "correlation_heatmap_mouse_all_samples")
    # save(plot_6x6_sn(rpm_mus_norm, "Mouse"), "correlation_heatmap_mouse_6x6_sn")
    # save(plot_3x6_lys_sn(rpm_mus_norm, "Mouse"), "correlation_heatmap_mouse_3x6_lys_sn")

    print("\nDone. Manuscript-faithful Suppl16 plots in", OUT_DIR)
