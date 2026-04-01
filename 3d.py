"""
Figure 3d — GSEA pathway enrichment heatmap (IFN-gamma)
Notebook-faithful reproduction of Figure3/GSEA_plot_enriched_pathways_two_columns.ipynb.
Outputs to: refactoring_roadmap/3d_plots/
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "3d_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = ROOT / "3d_csv"
CSV_DIR.mkdir(exist_ok=True)
NTVE_GSEA = ROOT / "resources" / "gsea_ifn-gamma" / "ntve" / "gsea_report_for_na_pos_1727164795980.tsv"
LYSATE_GSEA = ROOT / "resources" / "gsea_ifn-gamma" / "lysate" / "gsea_report_for_na_pos_1727164776729.tsv"

TOP_N = 15
TOP_N_PANEL = 10


if not NTVE_GSEA.exists() or not LYSATE_GSEA.exists():
    raise FileNotFoundError("Missing precomputed IFN-gamma GSEA reports under resources/gsea_ifn-gamma.")

ntve = pd.read_csv(NTVE_GSEA, sep="\t")
lysate = pd.read_csv(LYSATE_GSEA, sep="\t")

ntve_small = ntve.set_index("NAME")[["NES", "FDR q-val"]].sort_values("FDR q-val").head(TOP_N)
lysate_small = lysate.set_index("NAME")[["NES", "FDR q-val"]].sort_values("FDR q-val").head(TOP_N)
merged = (
    ntve_small.merge(
        lysate_small,
        on="NAME",
        how="outer",
        suffixes=("_NTVE", "_lysate"),
    )
    .sort_values("NES_lysate", ascending=False)
    .reset_index()
)
merged["NAME"] = [name.replace("HALLMARK_", "") for name in merged["NAME"]]


sorted_df = (
    merged.dropna(subset=["FDR q-val_lysate", "FDR q-val_NTVE"])
    .sort_values("FDR q-val_lysate", ascending=True, na_position="last")
    .head(TOP_N_PANEL)
    .copy()
)
sorted_df["lysate_rank"] = sorted_df["FDR q-val_lysate"].rank(method="min", na_option="keep", ascending=False)
sorted_df["SN_rank"] = sorted_df["FDR q-val_NTVE"].rank(method="min", na_option="keep", ascending=False)
sorted_df["neg_log_10_FDR_NTVE"] = -np.log10(sorted_df["FDR q-val_NTVE"] + 1e-3)
sorted_df["neg_log_10_FDR_lysate"] = -np.log10(sorted_df["FDR q-val_lysate"] + 1e-3)
sorted_df.to_csv(CSV_DIR / "3d_gsea_dotplot_data.csv", index=False)
plot_data = (
    sorted_df[["NAME", "FDR q-val_lysate", "FDR q-val_NTVE"]]
    .set_index("NAME")
    .rename(columns={"FDR q-val_lysate": "Lysate", "FDR q-val_NTVE": "NTVE"})
)
mask = plot_data.isna()

plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(40 / 25.4, 100 / 25.4))
sns.heatmap(
    plot_data,
    cmap="Reds_r",
    vmin=0,
    vmax=0.1,
    center=0.05,
    annot=True,
    fmt=".3f",
    linewidths=0.5,
    linecolor="#d0d0d0",
    cbar=True,
    cbar_kws={"label": "q-value"},
    mask=mask,
    ax=ax,
)
ax.set_facecolor("white")
ax.set_title("Heatmap of FDR q-values")
ax.set_xlabel("Conditions")
ax.set_ylabel("Gene Sets")
ax.tick_params(axis="x", rotation=0)
ax.tick_params(axis="y", rotation=0)

fig.tight_layout()
for ext in ("svg", "png"):
    out = OUT_DIR / f"3d_gsea_dotplot.{ext}"
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")
plt.close(fig)
print("Done.")
