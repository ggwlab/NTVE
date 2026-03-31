"""
Figure 1g — Enrichment factor distribution
Standalone reproduction of Figure1/fig1_enrichment_factor_distribution.ipynb
Outputs to: refactoring_roadmap/1g_plots/
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from ntvetools import load_and_preprocess_featurecounts, load_gtf_df, calculate_enrichment_factors

# ── Paths ─────────────────────────────────────────────────────────────────────
HARMONIZED_COUNTS_FILE = ROOT / "resources/harmonized_harmonized_gene_counts_rv_stranded.txt"
GTF_FILE               = ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
OUT_DIR = Path(__file__).parent / "1g_plots"
OUT_DIR.mkdir(exist_ok=True)

# Three matched SN / Lysate replicate pairs
SAMPLE_PAIRS = [
    (['24L006479'], ['24L006460']),
    (['24L006480'], ['24L006461']),
    (['24L006481'], ['24L006462']),
]

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading GTF...")
gtf_dict = load_gtf_df(str(GTF_FILE))
gene_id_to_biotype = gtf_dict["gene_id_to_biotype"]
gene_id_is_mt      = gtf_dict["gene_id_is_mt"]

print("Loading counts...")
full_result = load_and_preprocess_featurecounts(str(HARMONIZED_COUNTS_FILE))
rpm_df = full_result["rpm_ensg_only"]
print(f"  {len(rpm_df)} genes")

# ── Enrichment per replicate pair ─────────────────────────────────────────────
log_enrichments_list = []

for sn_list, lys_list in SAMPLE_PAIRS:
    # Keep only genes with RPM >= 1 in the lysate replicate
    lys_col = lys_list[0]
    keep = rpm_df.index[rpm_df[lys_col] >= 1]
    filtered_result = {**full_result, "rpm_ensg_only": rpm_df.loc[keep]}

    enrich = calculate_enrichment_factors(
        gene_id_is_mt=gene_id_is_mt,
        gene_id_to_biotype=gene_id_to_biotype,
        count_df=filtered_result,
        sn_list=sn_list,
        lysate_list=lys_list,
        pseudocount=0,
        verbose=False,
    )
    ef = np.array(enrich["protein_coding_enrichment_factors"])
    log_ef = np.log10(ef)
    log_ef = log_ef[np.isfinite(log_ef)]
    log_enrichments_list.append(log_ef)
    print(f"  {sn_list[0]} vs {lys_list[0]}: {len(log_ef):,} protein-coding genes")

# ── Plot ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})

all_vals = np.concatenate(log_enrichments_list)
bins = np.linspace(all_vals.min() - 0.05, all_vals.max() + 0.05, 71)
centers = (bins[:-1] + bins[1:]) / 2

densities = np.array([
    np.histogram(le, bins=bins, density=True)[0]
    for le in log_enrichments_list
])
mean_d = densities.mean(axis=0)
se_d   = densities.std(axis=0, ddof=1) / np.sqrt(len(densities))

fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(centers, mean_d, color="#1f77b4", lw=1.5, label=f"Mean (n=3 replicates)")
ax.fill_between(centers, mean_d - 1.96 * se_d, mean_d + 1.96 * se_d,
                alpha=0.3, color="#1f77b4", label="95% CI")
ax.axvline(0, color="black", ls="--", lw=1, alpha=0.5)
ax.set_xlabel("Log₁₀ Enrichment Factor (SN / Lysate)", fontsize=8)
ax.set_ylabel("Density", fontsize=8)
ax.set_title("Protein-coding gene enrichment\n(NTVE SN vs Lysate)", fontsize=8)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()

for ext in ("svg", "png"):
    out = OUT_DIR / f"1g_enrichment_distribution.{ext}"
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
    print(f"Saved: {out}")

plt.close(fig)
print("Done.")
