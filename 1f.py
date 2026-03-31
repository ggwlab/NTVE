"""
Figure 1f — NTVE vs Lysate scatter
Standalone reproduction of Figure1/1f_scatter_refactored.ipynb
Outputs to: refactoring_roadmap/1f_plots/
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from scipy.stats import pearsonr, spearmanr
from adjustText import adjust_text
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "Figure1"))
sys.path.insert(0, str(Path(__file__).parent.parent))

from ntvetools import load_and_preprocess_featurecounts, load_gtf_df

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
HARMONIZED_COUNTS_FILE = ROOT / "resources/harmonized_harmonized_gene_counts_rv_stranded.txt"
GTF_FILE = ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
OUT_DIR = Path(__file__).parent / "1f_plots"
OUT_DIR.mkdir(exist_ok=True)

SN_SAMPLES     = ['24L006479', '24L006480', '24L006481']
LYSATE_SAMPLES = ['24L006460', '24L006461', '24L006462']

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading GTF...")
gtf_dict = load_gtf_df(str(GTF_FILE))
gene_id_to_name    = gtf_dict["gene_id_to_name"]
gene_id_to_biotype = gtf_dict["gene_id_to_biotype"]
gene_id_is_mt      = gtf_dict["gene_id_is_mt"]

print("Loading counts...")
result = load_and_preprocess_featurecounts(str(HARMONIZED_COUNTS_FILE))
rpm_df = result["rpm_ensg_only"]
print(f"  {len(rpm_df)} genes loaded")

# ── Build gene table ──────────────────────────────────────────────────────────
df = pd.DataFrame(index=rpm_df.index)
df["sn_avg"]  = rpm_df[SN_SAMPLES].mean(axis=1)
df["lys_avg"] = rpm_df[LYSATE_SAMPLES].mean(axis=1)
df["GeneName"]          = df.index.map(lambda x: gene_id_to_name.get(x, x))
df["is_mt"]             = df.index.map(lambda x: gene_id_is_mt.get(x, False))
df["is_protein_coding"] = df.index.map(lambda x: gene_id_to_biotype.get(x, "") == "protein_coding")
df["enrichment_factor"] = (df["sn_avg"] + 0.1) / (df["lys_avg"] + 0.1)

# ── Plot ──────────────────────────────────────────────────────────────────────
matplotlib.rcParams.update({"text.usetex": False, "svg.fonttype": "none"})
plt.rc("font", size=8)

def log10p(x):
    return np.log10(x + 1e-3)

mt_mask    = df["is_mt"]
pc_mask    = df["is_protein_coding"] & ~df["is_mt"]
other_mask = ~df["is_protein_coding"] & ~df["is_mt"]

non_mt = df[~mt_mask]
pearson_r,  _ = pearsonr( log10p(non_mt["lys_avg"]), log10p(non_mt["sn_avg"]))
spearman_r, _ = spearmanr(log10p(non_mt["lys_avg"]), log10p(non_mt["sn_avg"]))

size_in = 82.5612343297975 / 25.4
fig, ax = plt.subplots(figsize=(size_in, size_in))
ax.set_position([0.25, 0.25, 0.65, 0.65])

ax.scatter(log10p(df[pc_mask]["lys_avg"]),    log10p(df[pc_mask]["sn_avg"]),
           c="#1f77b4", s=3, alpha=0.6, label="Protein-coding", rasterized=True)
ax.scatter(log10p(df[other_mask]["lys_avg"]), log10p(df[other_mask]["sn_avg"]),
           c="orange",  s=2, alpha=0.5, label="Noncoding",      rasterized=True)
ax.scatter(log10p(df[mt_mask]["lys_avg"]),    log10p(df[mt_mask]["sn_avg"]),
           c="red",     s=5, alpha=0.8, label="MT-encoded")

ax.set_xlabel("RPM per gene — Lysate", fontsize=8)
ax.set_ylabel("RPM per gene — SN",     fontsize=8)
ax.set_aspect("equal")
ax.annotate(f"Pearson r: {pearson_r:.2f}\nSpearman ρ: {spearman_r:.2f}\n(non-MT genes)",
            xy=(0.05, 1.05), xycoords="axes fraction", fontsize=8)

# Label very low SN points
low_sn = df[df["enrichment_factor"] < 1e-3].copy()
texts = []
for _, row in low_sn.iterrows():
    t = ax.text(log10p(row["lys_avg"]), log10p(row["sn_avg"]), row["GeneName"],
                fontsize=4, color="black")
    t.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground="white"),
                        path_effects.Normal()])
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
except Exception as e:
    print(f"adjust_text warning: {e}")

nice = np.concatenate(([0, 0.01, 0.1], 10**np.arange(0, 5)))
ticks = log10p(nice)
ax.set_xticks(ticks); ax.set_yticks(ticks)
fmt = lambda v: "n.d." if v == 0 else f"$10^{{{int(np.log10(v))}}}$"
ax.set_xticklabels([fmt(v) for v in nice], fontsize=8)
ax.set_yticklabels([fmt(v) for v in nice], fontsize=8)
ax.legend(fontsize=8, loc="upper left")

for ext in ("svg", "png"):
    out = OUT_DIR / f"1f_scatter.{ext}"
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
    print(f"Saved: {out}")

plt.close(fig)
print("Done.")
