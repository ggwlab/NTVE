"""
Figure 3d — GSEA pathway enrichment dotplot (IFN-γ)
Standalone reproduction of Figure3/GSEA_plot_enriched_pathways_two_columns.ipynb
Outputs to: refactoring_roadmap/3d_plots/
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent / "3d_plots"
OUT_DIR.mkdir(exist_ok=True)
FIG3 = ROOT / "Figure3"

# ── Load GSEA outputs ─────────────────────────────────────────────────────────
# Positive enrichment for NTVE and Lysate
pos_ntve   = pd.read_csv(sorted((FIG3 / "IFN-gamma_SN.GseaPreranked.1727164795980").glob("gsea_report_for_na_pos*.tsv"))[0], sep='\t')
pos_lysate = pd.read_csv(sorted((FIG3 / "IFN-gamma_lysate.GseaPreranked.1727164776729").glob("gsea_report_for_na_pos*.tsv"))[0], sep='\t')

top_n = 15
ntve_top   = pos_ntve.sort_values("FDR q-val").head(top_n).set_index("NAME")[["NES","FDR q-val"]]
lysate_top = pos_lysate.sort_values("FDR q-val").head(top_n).set_index("NAME")[["NES","FDR q-val"]]

merged = ntve_top.merge(lysate_top, on="NAME", how="outer", suffixes=("_NTVE","_lysate")).reset_index()
merged["NAME"] = merged["NAME"].str.replace("HALLMARK_","")

# sort by lysate FDR, NaN last
merged = merged.sort_values("FDR q-val_lysate", na_position='last')
merged["neg_log10_FDR_NTVE"]   = -np.log10(merged["FDR q-val_NTVE"]   + 1e-3)
merged["neg_log10_FDR_lysate"] = -np.log10(merged["FDR q-val_lysate"] + 1e-3)

# ── Plot ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({'font.size': 8, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(4, 6))

# Lysate (open circles)
scatter_lys = ax.scatter(merged["NES_lysate"], merged["NAME"],
                         s=merged["neg_log10_FDR_lysate"].fillna(0) * 20 + 10,
                         facecolors='none', edgecolors='#1f77b4', linewidths=1.2,
                         label='Lysate', zorder=3)
# NTVE (filled circles)
scatter_ntve = ax.scatter(merged["NES_NTVE"], merged["NAME"],
                          s=merged["neg_log10_FDR_NTVE"].fillna(0) * 20 + 10,
                          c='grey', edgecolors='black', linewidths=0.5,
                          label='NTVE', zorder=3)

ax.axvline(0, color='grey', lw=0.5, ls='--')
ax.set_xlabel('NES', fontsize=8)
ax.set_title('IFN-γ GSEA (Hallmark, top 15)\nSize = -log10(FDR)', fontsize=8)
ax.legend(fontsize=7)
fig.tight_layout()

for ext in ('svg', 'png'):
    out = OUT_DIR / f"3d_gsea_dotplot.{ext}"
    fig.savefig(out, format=ext, bbox_inches='tight', dpi=150)
    print(f"Saved: {out}")
plt.close(fig)
print("Done.")
