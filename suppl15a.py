"""
Supplementary Figure 15a — STAR-based sample correlation heatmaps (IFN-γ experiment)
Standalone reproduction of Suppl15/STAR_Correlations_IFNg-simplified.ipynb

Outputs to: refactoring_roadmap/suppl15a_plots/
  suppl15a_hexbin_lysate_vs_ntve.{svg,png}   — hexbin scatter (avg Lysate vs avg NTVE)

Notebook helper outputs remain available in code, but are disabled by default:
  suppl15a_correlation_pearson.{svg,png}
  suppl15a_correlation_spearman.{svg,png}
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import pearsonr, spearmanr

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

from ntvetools import load_gtf_df

ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent / "suppl15a_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "suppl15a_csv"
CSV_DIR.mkdir(exist_ok=True)

matplotlib.rcParams["text.usetex"] = False
matplotlib.rcParams["svg.fonttype"] = "none"
plt.rc("font", size=8)

# Sample mapping: biological name → sample ID
LYSATE_SAMPLES = {
    "Lysate Rep. 1": "24L006460",
    "Lysate Rep. 2": "24L006461",
    "Lysate Rep. 3": "24L006462",
}
NTVE_SAMPLES = {
    "NTVE Rep. 1": "24L006479",
    "NTVE Rep. 2": "24L006480",
    "NTVE Rep. 3": "24L006481",
}
# Averages used for hexbin (IFN-γ condition)
LYSATE_IFN_SAMPLES = ["24L006469", "24L006470", "24L006471"]
NTVE_IFN_SAMPLES   = ["24L006488", "24L006489", "24L006490"]

HEATMAP_SIZE_IN = 82.5612343297975 / 25.4   # 82.56 mm → inches


# --- Load reference data ---
print("Loading GTF data...")
gtf_data = load_gtf_df(str(ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_name  = gtf_data["gene_id_to_name"]
gene_id_is_mt    = gtf_data["gene_id_is_mt"]

# --- Load featureCounts table ---
print("Loading featureCounts table...")
fc = pd.read_csv(
    ROOT / "resources" / "harmonized_harmonized_gene_counts_rv_stranded.txt",
    sep="\t", skiprows=1,
).set_index("Geneid")

# Strip path + keep only the sample ID (first underscore-delimited token of filename)
sample_cols = sorted(fc.columns[5:])
keep_cols = [c for c in sample_cols if "toTranscriptome" not in c]
rename_map = {c: c.split("/")[-1].split("_")[0] for c in keep_cols}
sample_ids = list(rename_map.values())

# Select only sample columns (drops Chr/Start/End/Strand/Length metadata columns)
counts = fc[keep_cols].rename(columns=rename_map)

# --- RPM normalisation, filter to ENSG ---
print("Computing RPM...")
counts_ensg = counts[counts.index.str.startswith("ENSG")]
rpm = counts_ensg.div(counts_ensg.sum(axis=0) * 1e-6)

# --- Annotations ---
rpm["GeneName"] = [gene_id_to_name.get(g, g) for g in rpm.index]
rpm["is_mt"]    = [gene_id_is_mt.get(g, False)  for g in rpm.index]

# --- Averages and enrichment factor (for hexbin) ---
rpm["lys_avg"] = rpm[LYSATE_IFN_SAMPLES].mean(axis=1)
rpm["ntve_avg"] = rpm[NTVE_IFN_SAMPLES].mean(axis=1)
rpm["enrichment_factor"] = (rpm["ntve_avg"] + 0.1) / (rpm["lys_avg"] + 0.1)
rpm.reset_index().rename(columns={"index": "GeneID"}).to_csv(
    CSV_DIR / "suppl15a_hexbin_gene_points.csv", index=False
)

# --- Log10 columns for correlation analysis ---
log10p = lambda x: np.log10(x + 1e-3)
for label, sid in {**LYSATE_SAMPLES, **NTVE_SAMPLES}.items():
    rpm[label] = rpm[sid].apply(log10p)

CORR_COLS = list(LYSATE_SAMPLES.keys()) + list(NTVE_SAMPLES.keys())


# --- Helpers ---

def save_plot(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=300)
        print(f"Saved: {out}")


# --- Plot 1: Hexbin scatter ---

def plot_hexbin(df: pd.DataFrame, col_x: str, col_y: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(HEATMAP_SIZE_IN, HEATMAP_SIZE_IN))
    ax.set_position([0.25, 0.25, 0.65, 0.65])

    x_all = log10p(df[col_x])
    y_all = log10p(df[col_y])
    pearson_r, _  = pearsonr(x_all, y_all)
    spearman_r, _ = spearmanr(x_all, y_all)

    mt = df[df["is_mt"]]
    x_mt = log10p(mt[col_x])
    y_mt = log10p(mt[col_y])

    extent = [x_all.min(), x_all.max(), y_all.min(), y_all.max()]

    hb = ax.hexbin(x_all, y_all, gridsize=50, extent=extent,
                   cmap="inferno", norm=mcolors.LogNorm(vmax=1000))

    mt_cmap = mcolors.LinearSegmentedColormap.from_list(
        "", [(0, 0, 0, 0), (0.5, 0.5, 0.5, 1)]
    )

    class SaturatedNorm(mcolors.Normalize):
        def __call__(self, value, clip=None):
            return np.where(np.asarray(value) > 0, 1.0, 0.0)

    ax.hexbin(x_mt, y_mt, gridsize=50, extent=extent,
              cmap=mt_cmap, norm=SaturatedNorm(), mincnt=1,
              label="MT encoded genes")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = plt.colorbar(hb, cax=cax)
    cbar.set_label("Individual genes per bin", fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    ax.set_xlabel("CPM per gene for lysate", fontsize=8)
    ax.set_ylabel("CPM per gene for NTVE", fontsize=8)
    ax.set_aspect("equal")
    ax.annotate(
        f"Pearson r: {pearson_r:.2f}\nSpearman ρ: {spearman_r:.2f}\n(all genes)",
        xy=(0.05, 1.05), xycoords="axes fraction", fontsize=8, color="black",
    )

    # Label low-enrichment outliers
    low = df.query("enrichment_factor < 1e-2").copy()
    low["x"] = low[col_x].apply(log10p)
    low["y"] = low[col_y].apply(log10p)
    texts = []
    for _, row in low.iterrows():
        t = ax.text(row["x"], row["y"], row["GeneName"], fontsize=4, color="black")
        t.set_path_effects([
            path_effects.Stroke(linewidth=0.5, foreground="white"),
            path_effects.Normal(),
        ])
        texts.append(t)
    if HAS_ADJUSTTEXT and texts:
        try:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
        except Exception as e:
            print(f"Warning: adjustText failed: {e}")

    nice = np.concatenate(([0, 0.01, 0.1], 10 ** np.arange(0, 5)))
    ticks = log10p(nice)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    def fmt(v):
        return "n.d." if v == 0 else f"$10^{{{int(np.log10(v))}}}$"

    ax.set_xticklabels([fmt(v) for v in nice], fontsize=8)
    ax.set_yticklabels([fmt(v) for v in nice], fontsize=8)
    ax.legend(fontsize=8, loc="upper left")

    return fig


print("Plotting hexbin scatter...")
fig = plot_hexbin(rpm, "lys_avg", "ntve_avg")
save_plot(fig, "suppl15a_hexbin_lysate_vs_ntve")
plt.close(fig)


# --- Plot 2 & 3: Correlation heatmaps ---

def plot_correlation_heatmap(corr_matrix: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(HEATMAP_SIZE_IN, HEATMAP_SIZE_IN))
    ax.set_position([0.25, 0.25, 0.65, 0.65])
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 8
    sns.heatmap(
        corr_matrix,
        cmap="inferno",
        vmin=0.85,
        vmax=1.0,
        annot=True,
        fmt=".2f",
        square=True,
        cbar_kws={"shrink": 0.5},
        annot_kws={"size": 8},
        ax=ax,
    )
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    return fig


# Notebook helper outputs intentionally disabled by default:
#
# for method in ("pearson", "spearman"):
#     print(f"Plotting {method} correlation heatmap...")
#     corr = rpm[CORR_COLS].corr(method)
#     fig = plot_correlation_heatmap(corr)
#     save_plot(fig, f"suppl15a_correlation_{method}")
#     plt.close(fig)

print("Done.")
