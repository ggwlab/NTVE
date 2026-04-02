"""
Figure 1k — 3' end coverage bias: PolyA Lysate vs PolyA NTVE
Standalone reproduction of Suppl5/plot_3prime_bias_refactored.ipynb
(section "Execution: Lysate vs NTVE (PolyA-selected)")
Outputs to: Figure1/1k_plots/
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "1k_plots"
OUT_DIR.mkdir(exist_ok=True)

COV_DATA = ROOT / "resources" / "coverage_arrays"
GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

LYSATE_FILES = [
    "coverage_arrays_24L006460.npz",
    "coverage_arrays_24L006461.npz",
    "coverage_arrays_24L006462.npz",
]
NTVE_FILES = [
    "coverage_arrays_24L006479.npz",
    "coverage_arrays_24L006480.npz",
    "coverage_arrays_24L006481.npz",
]

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
    coverage_counts = np.zeros(window_size)

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
            coverage_counts += 1
        else:
            aligned_matrix[idx, -profile_length:] = normalized
            coverage_counts[-profile_length:] += 1

    return np.nanmean(aligned_matrix, axis=0), coverage_counts


def combine_replicates(file_paths, transcript_id_to_strand, window_size=5000):
    replicate_means = []
    replicate_counts = []
    for path in file_paths:
        print(f"  Processing {path.name}...")
        with np.load(path) as data:
            mean, counts = process_coverage_profiles(data, transcript_id_to_strand, window_size)
            replicate_means.append(mean)
            replicate_counts.append(counts)

    individual_profiles = np.stack(replicate_means)
    return {
        "grand_mean": np.mean(individual_profiles, axis=0),
        "individual_profiles": individual_profiles,
        "mean_counts": np.mean(np.stack(replicate_counts), axis=0),
        "n_replicates": len(file_paths),
        "window_size": window_size,
    }


def sd(individual_profiles: np.ndarray) -> np.ndarray:
    return np.std(individual_profiles, axis=0, ddof=1)


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
        panels_below = heights[i + 1:] if i + 1 < n_plots else []
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


print("Loading GTF...")
gtf_df = load_gtf_df(str(GTF_PATH))["gtf_df"]
transcript_id_to_strand = gtf_df.query("feature == 'transcript'").set_index("transcript_id")["strand"]

conditions = {
    "Lysate": [COV_DATA / f for f in LYSATE_FILES],
    "NTVE": [COV_DATA / f for f in NTVE_FILES],
}
colors = {"Lysate": "#888888", "NTVE": "#76c7dc"}

results = {}
for name, files in conditions.items():
    print(f"Processing {name}...")
    results[name] = combine_replicates(files, transcript_id_to_strand)

fig, (ax1, ax2) = setup_axes(axes_width_mm=50, axes_height_mm=[20, 20], n_plots=2)

x = np.arange(-results["Lysate"]["window_size"] + 1, 1)
for name, result in results.items():
    error = sd(result["individual_profiles"])
    ax1.fill_between(x, result["grand_mean"] - error, result["grand_mean"] + error, alpha=0.3, color=colors[name],
                     label=f'{name} (n={result["n_replicates"]})')
    ax1.plot(x, result["grand_mean"], color=colors[name], linewidth=1)
    ax2.plot(x, result["mean_counts"], color=colors[name], linewidth=1)

ax1.set_ylabel("Avg.\nnorm. depth")
ax1.set_xticklabels([])
ax1.grid(True, alpha=0.3)
ax1.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

ax2.set_ylabel("Num.\ntranscripts")
ax2.set_xlabel("Pos. rel. to 3' end of Ensembl transcript")
ax2.grid(True, alpha=0.3)

save(fig, "1k_coverage_profiles")
plt.close(fig)
print("Done.")
