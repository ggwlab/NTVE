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

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

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
    panel_specs = [
        (
            r"NTVE$_{\mathrm{PABP}}$ cell line",
            [("Stable Lysate (Avg)", "Lysate"), ("Stable Exported (Avg)", "NTVE")],
        ),
        (
            "transient expression",
            [("DJJT-cellular", "Lysate"), ("DJJT-gag", "NTVE")],
        ),
    ]

    palette = {"Standard": "#d9d9d9", "Mitochondrial": "#ffae1a"}
    fig, axes = plt.subplots(2, 1, figsize=(3.1, 4.2), sharex=True)

    for ax, (panel_title, rows) in zip(axes, panel_specs):
        panel_df = combined[combined["sample"].isin([sample for sample, _ in rows])].copy()
        panel_df["group"] = panel_df["sample"].map(dict(rows))
        order = ["Lysate", "NTVE"]

        sns.violinplot(
            data=panel_df[panel_df["type"] == "Standard"],
            x="log_counts",
            y="group",
            order=order,
            orient="h",
            density_norm="width",
            bw_adjust=0.5,
            inner=None,
            cut=0,
            linewidth=1,
            color=palette["Standard"],
            ax=ax,
        )

        sns.stripplot(
            data=panel_df[panel_df["type"] == "Mitochondrial"],
            x="log_counts",
            y="group",
            order=order,
            orient="h",
            size=4,
            alpha=0.75,
            jitter=0.16,
            color=palette["Mitochondrial"],
            edgecolor="none",
            ax=ax,
        )

        ax.set_title(panel_title, loc="left", pad=2, fontweight="bold", fontsize=10)
        ax.set_ylabel("")
        ax.set_xlabel(r"Log$_{10}$ (normalized CPM)")
        ax.grid(False)
        ax.set_xlim(-3.7, 4.2)
        ax.set_xticks([-3.0, -2.0, 0.0, 2.0, 4.0])
        ax.set_xticklabels(["n.d.", "-2", "0", "2", "4"])
        ax.tick_params(axis="y", length=3, width=1)
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_linewidth(1)

    legend_handles = [
        plt.Line2D([0], [0], color="#666666", lw=4),
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor=palette["Mitochondrial"], markeredgecolor="none"),
    ]
    axes[0].legend(
        legend_handles,
        ["nuclear", "mitochondrial"],
        loc="upper right",
        frameon=False,
        ncol=2,
        handlelength=0.8,
        handletextpad=0.4,
        columnspacing=0.8,
        bbox_to_anchor=(1.0, 1.42),
    )

    fig.tight_layout(h_pad=1.0)
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
