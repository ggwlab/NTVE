"""
Supplementary Figure 5b,c
- 5b: detected-gene Venn diagram
- 5c: upper coverage-profile panel for total-RNA Lysate / NTVE / HIV-gag

Standalone reproduction of the plotted Supplementary Figure 5 panels from:
- Suppl5/compare_STAR_vs_Minimap.ipynb
- Suppl5/plot_3prime_bias_refactored.ipynb

Legacy exploratory outputs from the old script are intentionally omitted.
Outputs to: refactoring_roadmap/suppl5_plots/
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).parent.parent

from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "suppl5_plots"
OUT_DIR.mkdir(exist_ok=True)


def first_existing_dir(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


S5_DATA = first_existing_dir(ROOT / "resources" / "sup5_comparison_data", ROOT / "Suppl5" / "comparison_data")
COV_DATA = ROOT / "resources" / "coverage_arrays"
GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

STAR_LYS_COLS = ["24L006460", "24L006461", "24L006462"]
STAR_SN_COLS = ["24L006479", "24L006480", "24L006481"]
MINIMAP_LYS_COLS = ["24L006460", "24L006461", "24L006462"]
MINIMAP_SN_COLS = ["24L006479", "24L006480", "24L006481"]

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "svg.fonttype": "none",
        "figure.dpi": 300,
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=300)
        print(f"Saved: {out}")


def process_coverage_profiles(coverage_arrays, transcript_id_to_strand, window_size=5000):
    transcript_ids = list(coverage_arrays.files)
    aligned_matrix = np.full((len(transcript_ids), window_size), np.nan)

    for idx, transcript_id in enumerate(transcript_ids):
        profile = coverage_arrays[transcript_id]
        strand = transcript_id_to_strand.get(transcript_id, "+")
        if strand == "-":
            profile = profile[::-1]

        norm_factor = profile.sum() / len(profile)
        if norm_factor <= 0:
            continue

        normalized = (profile / norm_factor).clip(0, 100)
        profile_length = len(normalized)
        if profile_length >= window_size:
            aligned_matrix[idx, :] = normalized[-window_size:]
        else:
            aligned_matrix[idx, -profile_length:] = normalized

    return aligned_matrix


def combine_replicates(file_paths, transcript_id_to_strand, window_size=5000):
    individual_profiles = []
    for path in file_paths:
        print(f"  Processing {path.name}...")
        with np.load(path) as data:
            individual_profiles.append(process_coverage_profiles(data, transcript_id_to_strand, window_size))

    replicate_means = np.array([np.nanmean(profile_matrix, axis=0) for profile_matrix in individual_profiles])
    return {
        "grand_mean": np.nanmean(replicate_means, axis=0),
        "individual_profiles": replicate_means,
        "window_size": window_size,
    }


def ci95_half_width(individual_profiles: np.ndarray) -> np.ndarray:
    n = individual_profiles.shape[0]
    if n < 2:
        return np.zeros(individual_profiles.shape[1])
    se = np.std(individual_profiles, axis=0, ddof=1) / np.sqrt(n)
    return t.ppf(0.975, df=n - 1) * se


def setup_axes(axes_width_mm=50, axes_height_mm=20, n_plots=1):
    left_margin_mm = 15
    right_margin_mm = 30
    top_margin_mm = 10
    bottom_margin_mm = 10
    middle_margin_mm = 5

    if isinstance(axes_height_mm, (list, tuple)):
        heights = list(axes_height_mm)
    else:
        heights = [axes_height_mm] * n_plots

    fig_width_mm = left_margin_mm + axes_width_mm + right_margin_mm
    fig_height_mm = bottom_margin_mm + sum(heights) + middle_margin_mm * (n_plots - 1) + top_margin_mm
    fig = plt.figure(figsize=(fig_width_mm / 25.4, fig_height_mm / 25.4))

    axes = []
    for i in range(n_plots):
        panels_below = heights[i + 1 :] if i + 1 < n_plots else []
        b_mm = bottom_margin_mm + sum(panels_below) + middle_margin_mm * len(panels_below)
        ax = fig.add_axes(
            [
                left_margin_mm / fig_width_mm,
                b_mm / fig_height_mm,
                axes_width_mm / fig_width_mm,
                heights[i] / fig_height_mm,
            ]
        )
        axes.append(ax)

    return fig, axes


print("Reproducing Suppl 5b: detected-gene Venn diagram...")
star_tpm = pd.read_csv(S5_DATA / "STAR_protein_coding_tpm_normalized.csv", index_col=0)[STAR_LYS_COLS + STAR_SN_COLS]
minimap_tpm = pd.read_csv(S5_DATA / "Minimap_protein_coding_tpm_normalized.csv", index_col=0)
common_genes = [gene for gene in minimap_tpm.index if gene in star_tpm.index]
star_tpm_common = star_tpm.loc[common_genes]
minimap_tpm_common = minimap_tpm.loc[common_genes]

star_lys_detected = set(star_tpm_common.index[(star_tpm_common[STAR_LYS_COLS] > 0).all(axis=1)])
star_sn_detected = set(star_tpm_common.index[(star_tpm_common[STAR_SN_COLS] > 0).all(axis=1)])
minimap_lys_detected = set(minimap_tpm_common.index[(minimap_tpm_common[MINIMAP_LYS_COLS] > 0).all(axis=1)])
minimap_sn_detected = set(minimap_tpm_common.index[(minimap_tpm_common[MINIMAP_SN_COLS] > 0).all(axis=1)])

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
venn2([minimap_lys_detected, minimap_sn_detected], set_labels=("Lysate", "NTVE"), ax=axes[0])
axes[0].set_title(
    "Minimap: Genes Detected in All Replicates\n(TPM > 0 in all 3 replicates)",
    fontsize=14,
    fontweight="bold",
)
venn2([star_lys_detected, star_sn_detected], set_labels=("Lysate", "NTVE"), ax=axes[1])
axes[1].set_title(
    "STAR: Genes Detected in All Replicates\n(TPM > 0 in all 3 replicates)",
    fontsize=14,
    fontweight="bold",
)
plt.tight_layout()
save(fig, "Venn_detected_genes")
plt.close(fig)


print("Reproducing Suppl 5c: coverage-profile upper panel...")
gtf_df = load_gtf_df(str(GTF_PATH))["gtf_df"]
transcript_id_to_strand = gtf_df.query("feature == 'transcript'").set_index("transcript_id")["strand"]

conditions = {
    "Lysate": ["coverage_arrays_23L010008.npz", "coverage_arrays_23L010009.npz", "coverage_arrays_23L010011.npz"],
    "NTVE": ["coverage_arrays_24L011708.npz", "coverage_arrays_24L011709.npz", "coverage_arrays_24L011710.npz"],
    "HIV-gag": ["coverage_arrays_24L011711.npz", "coverage_arrays_24L011712.npz", "coverage_arrays_24L011713.npz"],
}
colors = {"Lysate": "#888888", "NTVE": "#76c7dc", "HIV-gag": "#ffbc22"}
results = {
    name: combine_replicates([COV_DATA / f for f in files], transcript_id_to_strand)
    for name, files in conditions.items()
}

fig, (ax,) = setup_axes(axes_width_mm=50, axes_height_mm=[20], n_plots=1)
for name, result in results.items():
    x = np.arange(-result["window_size"] + 1, 1)
    mean = result["grand_mean"]
    ci = ci95_half_width(result["individual_profiles"])
    ax.plot(x, mean, color=colors[name], linewidth=1, label=name)
    ax.fill_between(x, mean - ci, mean + ci, color=colors[name], alpha=0.3)

ax.set_ylabel("Avg.\nnorm. depth")
ax.set_xlabel("Pos. rel. to 3' end of Ensembl transcript")
ax.grid(True, alpha=0.3)
ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
save(fig, "comparative_profiles_three_conditions_average_lysate")
plt.close(fig)

print("Done.")
