"""
Supplementary Figure 11a — Benchmarking comparison enrichment distributions
Standalone reproduction of Suppl11/enrichment_factor_distribution.ipynb
Outputs to: refactoring_roadmap/suppl11a_plots/
"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir

OUT_DIR = Path(__file__).parent / "suppl11a_plots"
OUT_DIR.mkdir(exist_ok=True)

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 8

NAME_MAPPING = {"S1": "NTVE", "S2": "HIV-gag WT", "S3": "TRACE-seq", "S4": "MMLV-gag WT", "S5": "HIV-gag:MCP"}


def calculate_enrichment_from_tpm(tpm_df: pd.DataFrame, sn_samples: list[str], lysate_samples: list[str], lysate_tpm_cutoff: float = 1.0) -> list[np.ndarray]:
    enrichment_factors_list = []
    for sn_sample, lysate_sample in zip(sn_samples, lysate_samples):
        sn_tpm = tpm_df[sn_sample].values
        lysate_tpm = tpm_df[lysate_sample].values
        mask = lysate_tpm >= lysate_tpm_cutoff
        enrichment = sn_tpm[mask] / lysate_tpm[mask]
        log_enrichment = np.log10(enrichment)
        enrichment_factors_list.append(log_enrichment[np.isfinite(log_enrichment)])
    return enrichment_factors_list


def plot_enrichment_histograms_tpm(enrichment_list: list[np.ndarray], title: str, bins: np.ndarray, xlim: tuple[float, float], ylim: tuple[float, float]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(4, 3))
    densities = [np.histogram(vals, bins=bins, density=True)[0] for vals in enrichment_list if len(vals) > 0]
    densities = np.array(densities)
    centers = (bins[:-1] + bins[1:]) / 2
    mean_density = np.mean(densities, axis=0)
    se_density = np.std(densities, axis=0, ddof=1) / np.sqrt(len(densities))
    ax.plot(centers, mean_density, color="#1f77b4", linewidth=1.5, label=f"Mean (n={len(densities)})")
    ax.fill_between(centers, mean_density - 1.96 * se_density, mean_density + 1.96 * se_density, alpha=0.3, color="#1f77b4", label="95% CI")
    ax.set_xlabel("Log10 Enrichment Factor (SN/Lysate)")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.axvline(x=0, color="black", linestyle="--", alpha=0.5, linewidth=1)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.tight_layout()
    return fig


def save_plot(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
        print(f"Saved: {out}")


print("Loading benchmarking comparison gene-level TPM matrix...")
tpm_df = pd.read_csv(get_generated_figure2_loaded_dir() / "elo_tpm_matrix_protein_coding.csv", index_col=0)
tpm_data = tpm_df.drop(columns=["GeneName"])
all_construct_data = {}
for construct in NAME_MAPPING:
    sn_cols = sorted([c for c in tpm_data.columns if c.startswith(f"SN_{construct}_")])
    lys_cols = sorted([c for c in tpm_data.columns if c.startswith(f"Lysate_{construct}_")])
    n_reps = min(len(sn_cols), len(lys_cols))
    all_construct_data[construct] = calculate_enrichment_from_tpm(tpm_data, sn_cols[:n_reps], lys_cols[:n_reps])
all_values = np.concatenate([np.concatenate(vals) for vals in all_construct_data.values()])
global_bins = np.linspace(np.min(all_values), np.max(all_values), 71)
global_y_max = 0.0
for enrichments in all_construct_data.values():
    densities = [np.histogram(vals, bins=global_bins, density=True)[0] for vals in enrichments if len(vals) > 0]
    densities = np.array(densities)
    ci_upper = np.mean(densities, axis=0) + 1.96 * np.std(densities, axis=0, ddof=1) / np.sqrt(len(densities))
    global_y_max = max(global_y_max, float(np.max(ci_upper)))
for construct, label in NAME_MAPPING.items():
    fig = plot_enrichment_histograms_tpm(all_construct_data[construct], title=f"Benchmarking comparison: {label}", bins=global_bins, xlim=(global_bins[0], global_bins[-1]), ylim=(0, global_y_max * 1.05))
    save_plot(fig, f"suppl11a_{construct.lower()}_enrichment_distribution")
    plt.close(fig)
print("Done.")
