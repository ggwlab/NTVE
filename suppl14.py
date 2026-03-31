"""
Supplementary Figure 14a-c — CRISPRa ASCL1 scatter and volcano plots.

Source notebook: Suppl14/Load_DE_From_R.ipynb

Outputs to: refactoring_roadmap/suppl14_plots/
  suppl14a_ascl1_filtered_only_significant
  suppl14_volcano_SN
  suppl14_volcano_lysate
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns
from adjustText import adjust_text
from pathlib import Path
from scipy.stats import pearsonr
from matplotlib.ticker import MaxNLocator
import sys

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "suppl14_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "suppl14_csv"
CSV_DIR.mkdir(exist_ok=True)
S14 = ROOT / "resources" / "suppl14"

gtf = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_name = gtf["gene_id_to_name"]
gene_id_to_gene_biotype = gtf["gene_id_to_biotype"]


def load_deseq(fname: str, suffix: str) -> pd.DataFrame:
    csv_path = S14 / fname
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Missing Suppl14 DESeq2 input: {csv_path}. "
            "Expected canonical static inputs under resources/suppl14/."
        )
    df = pd.read_csv(csv_path)
    df.rename(columns={df.columns[0]: "GeneID", "Geneid": "GeneID"}, inplace=True)
    keep = ["GeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()
    return df.rename(columns={c: f"{c}_{suffix}" for c in keep if c != "GeneID"})


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
        print(f"Saved: {out}")
    plt.close(fig)


def volcano_plot(
    merged_df: pd.DataFrame,
    log2fc_col: str,
    padj_col: str,
    gene_name_col: str,
    title: str,
    save_stem: str,
    genes_to_plot: list[int],
    p_cutoff: float = 1e-7,
) -> None:
    fig = plt.figure(figsize=(5, 5))
    axes_width_inch = 50 / 25.4
    axes_height_inch = 50 / 25.4
    ax = fig.add_axes([0.5 - axes_width_inch / 2 / 5, 0.5 - axes_height_inch / 2 / 5, axes_width_inch / 5, axes_height_inch / 5])

    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rc("font", size=8)

    df = merged_df.copy()
    df[padj_col] = df[padj_col] + 1e-300
    df["-log10(padj)"] = -np.log10(df[padj_col])
    df["Significant"] = df[padj_col] < p_cutoff

    sig_df = df[df["Significant"]]
    nonsig_df = df[~df["Significant"]]

    ax.scatter(
        nonsig_df[log2fc_col],
        nonsig_df["-log10(padj)"],
        color="grey",
        edgecolor="white",
        linewidth=0.5,
        alpha=1.0,
        s=15,
        rasterized=True,
    )
    ax.scatter(
        sig_df[log2fc_col],
        sig_df["-log10(padj)"],
        color="red",
        edgecolor="white",
        linewidth=0.5,
        alpha=1.0,
        s=15,
        rasterized=True,
    )

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", label=r"$\mathit{p}\geq10^{-7}$", markersize=6),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="red", label=r"$\mathit{p}<10^{-7}$", markersize=6),
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

    for x in [-5, -2.5, 0, 2.5, 5, 7.5]:
        alpha = 0.5 if x == 0 else 0.1
        ax.axvline(x=x, color="grey", linestyle="--", alpha=alpha)
    ax.axhline(y=-np.log10(p_cutoff), color="grey", linestyle="--", alpha=0.5)
    ax.set_xticks([-5, -2.5, 0, 2.5, 5, 7.5])
    ax.set_xlabel(r"Log$_{2}$ Fold Change")
    ax.set_ylabel(r"-Log$_{10}$(Adjusted $\mathit{p}$-value + $10^{-300}$)")
    ax.set_title(title, fontsize=10)

    top = df.loc[genes_to_plot]
    texts = []
    for _, row in top.iterrows():
        if pd.notna(row[gene_name_col]):
            text = ax.text(row[log2fc_col], -np.log10(row[padj_col]), row[gene_name_col], fontsize=6, color="black")
            text.set_path_effects([path_effects.Stroke(linewidth=1, foreground="white"), path_effects.Normal()])
            texts.append(text)
    if texts:
        try:
            adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
        except Exception as exc:
            print(f"adjust_text: {exc}")

    save(fig, save_stem)


def log_fc_scatter(
    merged_df: pd.DataFrame,
    save_stem: str,
    counting_level: str = "genes",
    filtering: str = "None",
    fc_annotation_cutoff_value: float = 5,
    axis_range=None,
) -> None:
    fig = plt.figure(figsize=(5, 5))
    axes_width_inch = 50 / 25.4
    axes_height_inch = 50 / 25.4
    ax = fig.add_axes([0.5 - axes_width_inch / 2 / 5, 0.5 - axes_height_inch / 2 / 5, axes_width_inch / 5, axes_height_inch / 5])

    plt.rcParams["svg.fonttype"] = "none"
    plt.rc("font", size=8)

    filtered_data = merged_df.dropna(
        subset=["log2FoldChange_lysate", "log2FoldChange_SN", "log_baseMean_lysate", "lysate_threshold", "SN_threshold"]
    ).copy()
    filtered_data["Log2 FoldChange in Lysate"] = filtered_data["minuslog2FoldChange_lysate"]
    filtered_data["Log2 FoldChange in NTVE"] = filtered_data["minuslog2FoldChange_SN"]

    filtered_data["threshold"] = "None"
    filtered_data.loc[(filtered_data["lysate_threshold"]) & (filtered_data["SN_threshold"]), "threshold"] = "Both"
    filtered_data.loc[(filtered_data["lysate_threshold"]) & (~filtered_data["SN_threshold"]), "threshold"] = "Lysate"
    filtered_data.loc[(~filtered_data["lysate_threshold"]) & (filtered_data["SN_threshold"]), "threshold"] = "NTVE"
    filtered_data = filtered_data.sort_values("threshold", ascending=False)

    color_dict = {"Both": "#8BC46F", "Lysate": "#FFE222", "NTVE": "#4EC1FF", "None": "white"}
    sns.scatterplot(
        data=filtered_data,
        x="Log2 FoldChange in Lysate",
        y="Log2 FoldChange in NTVE",
        hue="threshold",
        palette=color_dict,
        edgecolor="black",
        linewidth=0.5,
        s=25,
        ax=ax,
    )

    r_value_all, p_value_all = pearsonr(filtered_data["Log2 FoldChange in Lysate"], filtered_data["Log2 FoldChange in NTVE"])
    significant_data = filtered_data[filtered_data["threshold"] != "None"]
    r_value_sig, p_value_sig = pearsonr(significant_data["Log2 FoldChange in Lysate"], significant_data["Log2 FoldChange in NTVE"])

    high_de_points = filtered_data.query(
        "abs(log2FoldChange_lysate) > @fc_annotation_cutoff_value or abs(log2FoldChange_SN) > @fc_annotation_cutoff_value"
    ).copy()
    high_de_points["x"] = high_de_points["Log2 FoldChange in Lysate"]
    high_de_points["y"] = high_de_points["Log2 FoldChange in NTVE"]

    ax.set_aspect("equal")
    if axis_range is not None:
        min_val, max_val = axis_range
    else:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        min_val = min(xlim[0], ylim[0])
        max_val = max(xlim[1], ylim[1])
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.xaxis.set_major_locator(MaxNLocator(8, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(8, integer=True))
    ticks = sorted(set(ax.get_xticks()) | set(ax.get_yticks()))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.plot([min_val, max_val], [min_val, max_val], "--", color="grey", linewidth=1.5)
    ax.text(max_val - 0.05 * (max_val - min_val), max_val, "x=y", fontsize=8, color="grey")

    texts = []
    for _, row in high_de_points.iterrows():
        gene_name = row["gene_name_lysate"]
        if str(gene_name) != "nan":
            text = ax.text(row["x"], row["y"], gene_name, fontsize=6, color="black")
            text.set_path_effects([path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()])
            texts.append(text)
    if texts:
        try:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
        except Exception as exc:
            print(f"adjust_text: {exc}")

    textstr = "\n".join(
        (
            "All points:",
            f"Pearson r: {r_value_all:.2f}",
            f"p-value: {p_value_all:.2e}",
            "Significant points:",
            f"Pearson r: {r_value_sig:.2f}",
            f"p-value: {p_value_sig:.2e}",
            f"counting_level: {counting_level}",
            f"filtering: {filtering}",
        )
    )
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(1.1, 0.95, textstr, transform=ax.transAxes, fontsize=6, verticalalignment="top", bbox=props)
    ax.legend(title="Threshold", loc="center left", bbox_to_anchor=(1, 0.5))
    save(fig, save_stem)

    print(f"Pearson correlation coefficient (r): {r_value_all:.2f}")
    print(f"Pearson correlation p-value: {p_value_all:.2e}")
    print(f"Pearson correlation coefficient significant (r): {r_value_sig:.2f}")
    print(f"Pearson correlation p-value significant: {p_value_sig:.2e}")


def build_merged_df() -> pd.DataFrame:
    sn_df = load_deseq("SN_Dox_vs_Dox_ASCL1_base_mean_gt_0.csv", "SN")
    lys_df = load_deseq("Lysate_Dox_vs_Dox_ASCL1_base_mean_gt_0.csv", "lysate")

    merged_df = lys_df.merge(sn_df, on="GeneID", how="inner")
    merged_df["log_baseMean_lysate"] = np.log10(merged_df["baseMean_lysate"].replace(0, np.nan))
    merged_df["log_baseMean_SN"] = np.log10(merged_df["baseMean_SN"].replace(0, np.nan))
    merged_df["gene_name_lysate"] = [gene_id_to_name.get(gid, gid) for gid in merged_df["GeneID"]]
    merged_df["gene_name_SN"] = merged_df["gene_name_lysate"]
    merged_df["gene_biotype"] = [gene_id_to_gene_biotype.get(gid, "NA") for gid in merged_df["GeneID"]]
    merged_df["baseMean_SN_normalized"] = merged_df["baseMean_SN"] * 1e6 / merged_df["baseMean_SN"].sum()
    merged_df["baseMean_lysate_normalized"] = merged_df["baseMean_lysate"] * 1e6 / merged_df["baseMean_lysate"].sum()
    merged_df["minuslog2FoldChange_SN"] = -merged_df["log2FoldChange_SN"]
    merged_df["minuslog2FoldChange_lysate"] = -merged_df["log2FoldChange_lysate"]
    merged_df["SN_threshold"] = merged_df["padj_SN"] < 1e-7
    merged_df["lysate_threshold"] = merged_df["padj_lysate"] < 1e-7

    return merged_df


def main() -> None:
    merged_df = build_merged_df()
    merged_df.to_csv(CSV_DIR / "suppl14_ascl1_scatter_volcano_data.csv", index=False)
    genes_to_plot = list(
        merged_df.query(
            "padj_lysate < 1E-100 or padj_SN < 1E-150 or abs(log2FoldChange_SN) > 6.5 or (abs(log2FoldChange_lysate) > 6 and padj_lysate < 1E-7)"
        ).index
    )
    log_fc_scatter(
        merged_df.query("lysate_threshold == True or SN_threshold == True").copy(),
        save_stem="suppl14a_ascl1_filtered_only_significant",
        fc_annotation_cutoff_value=5,
    )
    volcano_plot(
        merged_df,
        "minuslog2FoldChange_SN",
        "padj_SN",
        "gene_name_SN",
        title="ASCL1 Transactivation: SN",
        save_stem="suppl14_volcano_SN",
        genes_to_plot=genes_to_plot,
    )
    volcano_plot(
        merged_df,
        "minuslog2FoldChange_lysate",
        "padj_lysate",
        "gene_name_SN",
        title="ASCL1 Transactivation: Lysate",
        save_stem="suppl14_volcano_lysate",
        genes_to_plot=genes_to_plot,
    )
    print("Done.")


if __name__ == "__main__":
    main()
