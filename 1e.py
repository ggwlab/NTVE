"""
Figure 1e — Mitochondrial vs nuclear encoded reads across sample types
Standalone reproduction of Figure1/Colab_MT_read_count_proportion_refactored.ipynb
Outputs to: refactoring_roadmap/1e_plots/
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "resources" / "ntvetools"))

from ntvetools import load_gtf_df

GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
COUNTS_PATH = ROOT / "resources" / "harmonized_harmonized_gene_counts_rv_stranded.txt"
OUT_DIR = Path(__file__).parent / "1e_plots"
OUT_DIR.mkdir(exist_ok=True)

matplotlib.rcParams["text.usetex"] = False
matplotlib.rcParams["svg.fonttype"] = "none"
plt.rc("font", size=10)


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=300)
        print(f"Saved: {out}")


def create_mt_violin_plot(non_mt_df: pd.DataFrame, mt_df: pd.DataFrame, sample_map: dict[str, str]) -> plt.Figure:
    def melt_data(df: pd.DataFrame, label: str) -> pd.DataFrame:
        melted = df.reset_index().melt(
            id_vars=["GeneID"],
            value_vars=list(sample_map.keys()),
            value_name="read_counts",
            var_name="sample",
        )
        melted["type"] = label
        return melted

    combined = pd.concat(
        [
            melt_data(non_mt_df, "Standard"),
            melt_data(mt_df, "Mitochondrial"),
        ],
        ignore_index=True,
    )
    combined["log_counts"] = np.log10(combined["read_counts"] + 1e-3)
    combined["sample_label"] = combined["sample"].map(sample_map)

    fig, ax = plt.subplots(figsize=(12, 6))
    palette = {"Standard": "lightgray", "Mitochondrial": "orange"}

    sns.violinplot(
        x="sample_label",
        y="log_counts",
        hue="type",
        data=combined[combined["type"] == "Standard"],
        hue_order=["Standard", "Mitochondrial"],
        density_norm="width",
        bw_adjust=0.5,
        alpha=0.3,
        inner="box",
        order=list(sample_map.values()),
        palette=palette,
        dodge=True,
        ax=ax,
    )

    sns.stripplot(
        x="sample_label",
        y="log_counts",
        hue="type",
        data=combined[combined["type"] == "Mitochondrial"],
        hue_order=["Standard", "Mitochondrial"],
        palette=palette,
        size=5,
        alpha=0.7,
        jitter=0.2,
        order=list(sample_map.values()),
        dodge=True,
        ax=ax,
    )

    ax.set_title("Mitochondrial vs Nuclear Encoded Reads\nAcross Sample Types", pad=20)
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) >= 4:
        ax.legend(
            [handles[0], handles[3]],
            ["Nuclear (Violin)", "Mitochondrial (Dots)"],
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
        )
    else:
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    ax.set_ylabel("Log10(Read Counts)")
    ax.set_xlabel("")
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha("right")
    tick_colors = ["darkblue", "darkred"] * 2
    for tick, color in zip(ax.get_xticklabels(), tick_colors):
        tick.set_color(color)

    fig.tight_layout()
    return fig


print("Loading GTF...")
gtf_data = load_gtf_df(str(GTF_PATH))
gtf_df = gtf_data["gtf_df"]

print("Loading raw counts...")
raw_counts_df = pd.read_csv(COUNTS_PATH, sep="\t", skiprows=1).rename(columns={"Geneid": "GeneID"})

print("Preparing RPM matrix...")
gene_metadata = gtf_df.drop_duplicates("gene_id").set_index("gene_id")
gene_id_to_is_mito = (gene_metadata["seqname"] == "hs_MT").to_dict()

sample_columns = raw_counts_df.columns[6:]
column_rename_map = {
    col: col.replace("./star_alignments/Project_NTVE/star_harmonized/", "").replace(
        "_untrimmed_harmonizedAligned.sortedByCoord.out.bam", ""
    )
    for col in sample_columns
}
processed_counts = raw_counts_df.rename(columns=column_rename_map).set_index("GeneID")
clean_sample_cols = [column_rename_map[c] for c in sample_columns]

rpm_df = processed_counts[clean_sample_cols].copy()
rpm_df = rpm_df / rpm_df.sum(axis=0) * 1e6
human_rpm_df = rpm_df.loc[rpm_df.index.str.startswith("ENSG")].copy()

human_rpm_df["Stable Lysate (Avg)"] = human_rpm_df[["24L006463", "24L006464", "24L006465"]].mean(axis=1)
human_rpm_df["Stable Exported (Avg)"] = human_rpm_df[["24L006482", "24L006483", "24L006484"]].mean(axis=1)

plot_cols_map = {
    "DJJT-cellular": "Transient Lysate",
    "DJJT-gag": "Transient Exported",
    "Stable Lysate (Avg)": "Stable Lysate",
    "Stable Exported (Avg)": "Stable Exported (NTVE)",
}

is_mito_mask = human_rpm_df.index.map(lambda x: gene_id_to_is_mito.get(x, False))
mt_df = human_rpm_df[is_mito_mask].copy()
non_mt_df = human_rpm_df[~is_mito_mask].copy()

print("Plotting mitochondrial vs nuclear read distributions...")
fig = create_mt_violin_plot(non_mt_df, mt_df, plot_cols_map)
save(fig, "MT_Content_Across_samples_new")
plt.close(fig)
print("Done.")
