"""
Figure 5b — Species-specific transcriptome recovery pie charts
Standalone reproduction of Figure5/Revisited_piecharts.ipynb
Outputs to: refactoring_roadmap/5b_plots/
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent / "5b_plots"
OUT_DIR.mkdir(exist_ok=True)
COUNTS_FILE = ROOT / "resources/harmonized_harmonized_gene_counts_rv_stranded.txt"

MAPPING = {
    "24L010907": "N2a-1",
    "24L010908": "N2a-2",
    "24L010909": "N2a-3",
    "24L010910": "HEK-1",
    "24L010911": "HEK-2",
    "24L010912": "HEK-3",
}

print("Loading counts...")
featurecounts = pd.read_csv(COUNTS_FILE, sep="\t", skiprows=1).set_index("Geneid")
sample_columns = sorted(featurecounts.columns[5:])
raw_counts_df = featurecounts.rename(columns={col: col.split("/")[-1].split("_")[0] for col in sample_columns})

gene_level_df = raw_counts_df.loc[raw_counts_df.index.str.contains("MUS") | raw_counts_df.index.str.contains("ENSG")]
gene_level_df_rpm = gene_level_df[list(MAPPING.keys())] / (gene_level_df[list(MAPPING.keys())].sum(axis=0) * 1e-6)
gene_level_df_rpm_mus = gene_level_df_rpm.loc[gene_level_df_rpm.index.str.contains("MUS")]
gene_level_df_rpm_homo = gene_level_df_rpm.loc[gene_level_df_rpm.index.str.contains("ENSG")]

hek_samples = ["24L010910", "24L010911", "24L010912"]
n2a_samples = ["24L010907", "24L010908", "24L010909"]
hek_data = [[gene_level_df_rpm_homo[sample].sum(), gene_level_df_rpm_mus[sample].sum()] for sample in hek_samples]
n2a_data = [[gene_level_df_rpm_homo[sample].sum(), gene_level_df_rpm_mus[sample].sum()] for sample in n2a_samples]

plt.rcParams.update({"svg.fonttype": "none", "font.family": "Arial", "font.size": 8})
pie_diameter_inches = 23 / 25.4
labels = ["Human", "Mouse"]
colors = ["#c94b45ff", "#ffd21eff"]
fig, axs = plt.subplots(2, 3, figsize=(pie_diameter_inches * 3.5, pie_diameter_inches * 2.5))

for i, data in enumerate(hek_data):
    axs[0, i].pie(data, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90)
    axs[0, i].set_title(f"HEK-{i+1}")

for i, data in enumerate(n2a_data):
    axs[1, i].pie(data, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90)
    axs[1, i].set_title(f"N2a-{i+1}")

fig.tight_layout()
for ext in ("svg", "png"):
    out = OUT_DIR / f"5b_pie_charts_species_distribution.{ext}"
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
    print(f"Saved: {out}")
plt.close(fig)
print("Done.")
