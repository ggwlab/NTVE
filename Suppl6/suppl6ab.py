"""
Supplementary Figure 6a, 6b — Export ratio distributions
(a) MRIi003-A hiPSCs
(b) MRIi003-A hiPSC-derived cardiomyocytes

Requires suppl6_plots/suppl6_cm_tpm_matrix.csv — run suppl6_data.py first.

The line shows the mean density of three replicates; the grey shading is the
95% confidence interval (as described in the figure legend).
Outputs to: refactoring_roadmap/suppl6_plots/
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = Path(__file__).parent / "suppl6_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "suppl6_csv"
CSV_DIR.mkdir(exist_ok=True)
MATRIX_FILE = OUT_DIR / "suppl6_cm_tpm_matrix.csv"

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 8

# Sample groups: (SN replicate, matched Lysate replicate)
PANELS = {
    "a": {
        "sn":     ["cm_ipsc_0_1", "cm_ipsc_0_2", "cm_ipsc_0_3"],
        "lysate": ["cm_ipsc_Lysate_1", "cm_ipsc_Lysate_2", "cm_ipsc_Lysate_3"],
        "title":  "hiPSC: Export ratio distribution\n(protein-coding, non-MT genes)",
        "stem":   "suppl6a_ipsc_export_ratio",
    },
    "b": {
        "sn":     ["cm_10_1", "cm_10_2", "cm_10_3"],
        "lysate": ["cm_cm_lysate_1", "cm_cm_lysate_2", "cm_cm_lysate_3"],
        "title":  "Cardiomyocyte: Export ratio distribution\n(protein-coding, non-MT genes)",
        "stem":   "suppl6b_cm_export_ratio",
    },
}

LYSATE_TPM_CUTOFF = 1.0   # exclude genes with lysate TPM < 1 (matches notebook)
N_BINS = 70


# --- Helper functions ---

def calculate_enrichment(
    tpm_df: pd.DataFrame,
    sn_samples: list[str],
    lysate_samples: list[str],
) -> list[np.ndarray]:
    """Return one log10(SN/Lysate) array per replicate pair, filtered to finite values."""
    enrichments = []
    for sn, lys in zip(sn_samples, lysate_samples):
        sn_tpm = tpm_df[sn].values
        lys_tpm = tpm_df[lys].values
        mask = lys_tpm >= LYSATE_TPM_CUTOFF
        log_enr = np.log10(sn_tpm[mask] / lys_tpm[mask])
        enrichments.append(log_enr[np.isfinite(log_enr)])
    return enrichments


def plot_enrichment(
    enrichments: list[np.ndarray],
    title: str,
    bins: np.ndarray,
    ylim: float,
) -> plt.Figure:
    densities = np.array(
        [np.histogram(e, bins=bins, density=True)[0] for e in enrichments if len(e) > 0]
    )
    centers = (bins[:-1] + bins[1:]) / 2
    mean_d = np.mean(densities, axis=0)
    se_d = np.std(densities, axis=0, ddof=1) / np.sqrt(len(densities))

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(centers, mean_d, color="#1f77b4", linewidth=1.5, label=f"Mean (n={len(densities)})")
    ax.fill_between(
        centers,
        mean_d - 1.96 * se_d,
        mean_d + 1.96 * se_d,
        alpha=0.3,
        color="#aec7e8",
        label="95% CI",
    )
    ax.axvline(x=0, color="black", linestyle="--", alpha=0.5, linewidth=1)
    ax.set_xlabel("Log₁₀ Export Ratio (SN / Lysate)")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=7)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, ylim)
    fig.tight_layout()
    return fig


def save_plot(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
        print(f"Saved: {out}")


# --- Main ---
if not MATRIX_FILE.exists():
    raise FileNotFoundError(
        f"{MATRIX_FILE} not found — run suppl6_data.py first."
    )

print("Loading TPM matrix...")
tpm_matrix = pd.read_csv(MATRIX_FILE, index_col=0)
tpm_data = tpm_matrix.drop(columns=["GeneName"])

# Pre-compute enrichments for both panels
all_enrichments = {
    panel: calculate_enrichment(tpm_data, cfg["sn"], cfg["lysate"])
    for panel, cfg in PANELS.items()
}

# Shared axis limits across panels (consistent comparison)
all_values = np.concatenate([e for enrs in all_enrichments.values() for e in enrs])
global_bins = np.linspace(
    np.percentile(all_values, 0.5),
    np.percentile(all_values, 99.5),
    N_BINS + 1,
)
global_ylim = 0.0
for enrs in all_enrichments.values():
    densities = np.array(
        [np.histogram(e, bins=global_bins, density=True)[0] for e in enrs if len(e) > 0]
    )
    ci_upper = np.mean(densities, axis=0) + 1.96 * np.std(densities, axis=0, ddof=1) / np.sqrt(len(densities))
    global_ylim = max(global_ylim, float(np.max(ci_upper)))

# Plot panels
raw_rows = []
curve_rows = []
for panel, cfg in PANELS.items():
    for rep_idx, enrichment_values in enumerate(all_enrichments[panel], start=1):
        raw_rows.extend(
            {
                "panel": panel,
                "stem": cfg["stem"],
                "replicate_pair": rep_idx,
                "log10_export_ratio": value,
            }
            for value in enrichment_values
        )
    densities = np.array(
        [np.histogram(e, bins=global_bins, density=True)[0] for e in all_enrichments[panel] if len(e) > 0]
    )
    centers = (global_bins[:-1] + global_bins[1:]) / 2
    mean_d = np.mean(densities, axis=0)
    se_d = np.std(densities, axis=0, ddof=1) / np.sqrt(len(densities))
    curve_rows.extend(
        {
            "panel": panel,
            "stem": cfg["stem"],
            "n_replicates": len(densities),
            "bin_center": center,
            "mean_density": mean,
            "se_density": se,
            "ci95_lower": mean - 1.96 * se,
            "ci95_upper": mean + 1.96 * se,
        }
        for center, mean, se in zip(centers, mean_d, se_d)
    )
    fig = plot_enrichment(
        all_enrichments[panel], cfg["title"], global_bins, global_ylim * 1.05
    )
    save_plot(fig, cfg["stem"])
    plt.close(fig)

pd.DataFrame(raw_rows).to_csv(CSV_DIR / "suppl6_export_ratio_raw_values.csv", index=False)
pd.DataFrame(curve_rows).to_csv(CSV_DIR / "suppl6_export_ratio_curve_summary.csv", index=False)
print("Done.")
