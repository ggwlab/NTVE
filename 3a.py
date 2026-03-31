"""
Figure 3a — IFN-gamma perturbation: filtered log fold-change scatter
Faithful standalone reproduction of the final `log_FC_scatter(...)` path used in
Figure3/Load_DE_From_R_IFNg.ipynb for:
merged_df.query("lysate_threshold == True or SN_threshold == True")
Outputs to: refactoring_roadmap/3a_plots/
"""

from pathlib import Path
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from scipy.stats import pearsonr
from scipy import stats

ROOT = Path(__file__).parent.parent


def first_existing_dir(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


FIG3 = first_existing_dir(ROOT / "resources" / "fig3", ROOT / "Figure3")
OUT_DIR = Path(__file__).parent / "3a_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "3a_csv"
CSV_DIR.mkdir(exist_ok=True)

sys.path.insert(0, str(ROOT))
from ntvetools import load_gtf_df

warnings.filterwarnings("ignore")
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rc("font", size=8)


def load_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df.rename(columns={"Geneid": "GeneID"})


def log_fc_scatter(
    merged_df: pd.DataFrame,
    title: str | None = None,
    save_title: str | None = None,
    counting_level: str = "genes",
    filtering: str = "None",
    fc_annotation_cutoff_value: float = 3,
    axis_range: tuple[float, float] | None = None,
) -> plt.Figure:
    axes_width_inch = 50 / 25.4
    axes_height_inch = 50 / 25.4
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([
        0.5 - axes_width_inch / 2 / 5,
        0.5 - axes_height_inch / 2 / 5,
        axes_width_inch / 5,
        axes_height_inch / 5,
    ])

    filtered_data = merged_df.dropna(
        subset=[
            "log2FoldChange_lysate",
            "log2FoldChange_SN",
            "log_baseMean_lysate",
            "lysate_threshold",
            "SN_threshold",
        ]
    ).copy()
    filtered_data["Log2 FoldChange in Lysate"] = -filtered_data["log2FoldChange_lysate"]
    filtered_data["Log2 FoldChange in NTVE"] = -filtered_data["log2FoldChange_SN"]

    filtered_data["threshold"] = "None"
    filtered_data.loc[
        (filtered_data["lysate_threshold"] == True) & (filtered_data["SN_threshold"] == True),
        "threshold",
    ] = "Both"
    filtered_data.loc[
        (filtered_data["lysate_threshold"] == True) & (filtered_data["SN_threshold"] != True),
        "threshold",
    ] = "Lysate"
    filtered_data.loc[
        (filtered_data["lysate_threshold"] != True) & (filtered_data["SN_threshold"] == True),
        "threshold",
    ] = "NTVE"
    filtered_data = filtered_data.sort_values("threshold", ascending=False)

    color_dict = {"Both": "#8BC46F", "Lysate": "#FFE222", "NTVE": "#4EC1FF", "None": "white"}

    print(f"Data points: {len(filtered_data)}")
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

    if title:
        plt.title(title)

    r_value_all, p_value_all = pearsonr(
        filtered_data["Log2 FoldChange in Lysate"],
        filtered_data["Log2 FoldChange in NTVE"],
    )

    significant_data = filtered_data[filtered_data["threshold"] != "None"]
    x = significant_data["Log2 FoldChange in Lysate"]
    y = significant_data["Log2 FoldChange in NTVE"]
    slope, intercept, r_value_sig, p_value_sig, std_err = stats.linregress(x, y)

    x_line = np.linspace(min(x), max(x), 100)
    y_line = slope * x_line + intercept

    n = len(x)
    mean_x = np.mean(x)
    t = stats.t.ppf(0.975, n - 2)
    mse = np.sum((y - (slope * x + intercept)) ** 2) / (n - 2)
    s_err_pred = np.sqrt(mse * (1 + 1 / n + (x_line - mean_x) ** 2 / np.sum((x - mean_x) ** 2)))
    y_upper = y_line + t * s_err_pred
    y_lower = y_line - t * s_err_pred

    ax.plot(x_line, y_line, color="red", linewidth=1.5, label="Regression line")
    ax.fill_between(x_line, y_lower, y_upper, color="gray", alpha=0.4, label="95% Prediction Interval")

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
        if str(row["gene_name_lysate"]) != "nan":
            text = ax.text(row["x"], row["y"], row["gene_name_lysate"], fontsize=6, color="black")
            text.set_path_effects([
                path_effects.Stroke(linewidth=2, foreground="white"),
                path_effects.Normal(),
            ])
            texts.append(text)
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

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
    fig.tight_layout()
    if save_title:
        fig.savefig(save_title, dpi=300, bbox_inches="tight")

    print(f"Pearson correlation coefficient (r): {r_value_all:.2f}")
    print(f"Pearson correlation p-value: {p_value_all:.2e}")
    print(f"Pearson correlation coefficient at least one below threhold (r): {r_value_sig:.2f}")
    print(f"Pearson correlation p-value at least one below threhold: {p_value_sig:.2e}")
    return fig


print("Loading DESeq2 results...")
lysate = load_results(FIG3 / "Dox_vs_IFNg_lysate_base_mean_gt_0.csv")
SN = load_results(FIG3 / "Dox_vs_IFNg_SN_base_mean_gt_0.csv")
merged_df = lysate.merge(SN, left_on="GeneID", right_on="GeneID", suffixes=["_lysate", "_SN"])

gtf_data = load_gtf_df(path_to_gtf=str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_name = gtf_data["gene_id_to_name"]
gene_id_to_gene_biotype = gtf_data["gene_id_to_biotype"]

merged_df["log_baseMean_lysate"] = np.log10(merged_df["baseMean_lysate"])
merged_df["log_baseMean_SN"] = np.log10(merged_df["baseMean_SN"])
merged_df["gene_name_lysate"] = [gene_id_to_name[gid] for gid in merged_df["GeneID"]]
merged_df["gene_name_SN"] = [gene_id_to_name[gid] for gid in merged_df["GeneID"]]
merged_df["gene_biotype"] = [gene_id_to_gene_biotype[gid] for gid in merged_df["GeneID"]]
merged_df["baseMean_SN_normalized"] = [val * 1e6 / merged_df["baseMean_SN"].sum() for val in merged_df["baseMean_SN"]]
merged_df["baseMean_lysate_normalized"] = [val * 1e6 / merged_df["baseMean_lysate"].sum() for val in merged_df["baseMean_lysate"]]
merged_df["SN_threshold"] = [True if val < 1e-7 else False for val in merged_df["padj_SN"]]
merged_df["lysate_threshold"] = [True if val < 1e-7 else False for val in merged_df["padj_lysate"]]

plot_df = merged_df.query("lysate_threshold == True or SN_threshold == True")
plot_df = plot_df.copy()
plot_df["plot_log2fc_lysate"] = -plot_df["log2FoldChange_lysate"]
plot_df["plot_log2fc_sn"] = -plot_df["log2FoldChange_SN"]
plot_df["plot_threshold_class"] = np.select(
    [
        plot_df["lysate_threshold"] & plot_df["SN_threshold"],
        plot_df["lysate_threshold"] & ~plot_df["SN_threshold"],
        ~plot_df["lysate_threshold"] & plot_df["SN_threshold"],
    ],
    ["Both", "Lysate", "NTVE"],
    default="None",
)
plot_df.to_csv(CSV_DIR / "3a_logFC_scatter_points.csv", index=False)
fig = log_fc_scatter(
    plot_df,
    fc_annotation_cutoff_value=4,
    save_title=str(OUT_DIR / "3a_logFC_scatter.svg"),
)
fig.savefig(OUT_DIR / "3a_logFC_scatter.png", dpi=300, bbox_inches="tight")
print(f"Saved: {OUT_DIR / '3a_logFC_scatter.png'}")
plt.close(fig)
print("Done.")
