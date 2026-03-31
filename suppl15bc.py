"""
Suppl Fig 15b,c — IFN-gamma volcano plots (Lysate and NTVE SN)
Standalone reproduction of Suppl15/Load_DE_From_R_IFNg.ipynb
Outputs to: refactoring_roadmap/suppl15bc_plots/
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from adjustText import adjust_text
from pathlib import Path
import sys

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "suppl15bc_plots"
OUT_DIR.mkdir(exist_ok=True)


def first_existing_dir(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


S15 = first_existing_dir(ROOT / "resources" / "fig3", ROOT / "Suppl15")

gtf = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_name = gtf["gene_id_to_name"]


def load_deseq(fname: str, suffix: str) -> pd.DataFrame:
    df = pd.read_csv(S15 / fname)
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
    df: pd.DataFrame,
    log2fc_col: str,
    padj_col: str,
    gene_name_col: str,
    save_stem: str,
    genes_to_plot: list[int],
    p_cutoff: float = 1e-7,
) -> None:
    fig = plt.figure(figsize=(5, 5))
    axes_width_inch = 50 / 25.4
    axes_height_inch = 50 / 25.4
    ax = fig.add_axes([0.5 - axes_width_inch / 2 / 5, 0.5 - axes_height_inch / 2 / 5, axes_width_inch / 5, axes_height_inch / 5])

    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rc("font", size=8)

    df = df.copy()
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

    print(f"Total genes: {len(df)}")
    print(f"Significant genes: {len(sig_df)}")
    print(f"Up-regulated: {len(sig_df[sig_df[log2fc_col] > 0])}")
    print(f"Down-regulated: {len(sig_df[sig_df[log2fc_col] < 0])}")


def build_merged_df() -> pd.DataFrame:
    lysate = load_deseq("Dox_vs_IFNg_lysate_base_mean_gt_0.csv", "lysate")
    sn = load_deseq("Dox_vs_IFNg_SN_base_mean_gt_0.csv", "SN")
    merged_df = lysate.merge(sn, on="GeneID", suffixes=["_lysate", "_SN"])
    merged_df["gene_name_SN"] = [gene_id_to_name.get(gid, gid) for gid in merged_df["GeneID"]]
    merged_df["minuslog2FoldChange_lysate"] = -merged_df["log2FoldChange_lysate"]
    merged_df["minuslog2FoldChange_SN"] = -merged_df["log2FoldChange_SN"]
    return merged_df


def main() -> None:
    merged_df = build_merged_df()
    genes_to_plot = list(
        merged_df.query("padj_lysate < 1E-150 or padj_SN < 1E-150 or abs(log2FoldChange_SN) > 6.5 or abs(log2FoldChange_lysate) > 6.5").index
    )
    if 47438 in genes_to_plot:
        genes_to_plot.remove(47438)

    volcano_plot(
        merged_df,
        "minuslog2FoldChange_lysate",
        "padj_lysate",
        "gene_name_SN",
        save_stem="suppl15b_volcano_lysate",
        genes_to_plot=genes_to_plot,
    )
    volcano_plot(
        merged_df,
        "minuslog2FoldChange_SN",
        "padj_SN",
        "gene_name_SN",
        save_stem="suppl15c_volcano_SN",
        genes_to_plot=genes_to_plot,
    )
    print("Done.")


if __name__ == "__main__":
    main()
