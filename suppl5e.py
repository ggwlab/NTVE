"""
Supplementary Figure 5e — Total-RNA all-human-transcript TPM by transcript length
Standalone reproduction of the all-human-transcript total-RNA plot block from
Figure1/length_resolved_analysis_refactored.ipynb.

Outputs to: refactoring_roadmap/suppl5e_plots/
"""

from pathlib import Path
import glob
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent

from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "suppl5e_plots"
OUT_DIR.mkdir(exist_ok=True)

GTF_FILE = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
QUANT_GLOB = str(ROOT / "resources" / "salmon_harmonized" / "*" / "quant.sf.gz")

TOTAL_RNA_COLS = [
    "24L011708", "24L011709", "24L011710",
    "24L011711", "24L011712", "24L011713",
    "23L010004", "23L010007", "23L010008",
    "23L010009", "23L010011", "23L010014",
    "24L006460", "24L006461", "24L006462",
]
TOTAL_RNA_GAG_PABP = ["24L011708", "24L011709", "24L011710"]
TOTAL_RNA_GAG_ONLY = ["24L011711", "24L011712", "24L011713"]
TOTAL_RNA_REFERENCE = ["23L010008", "23L010009", "23L010011"]

BINS = [0] + list(range(500, 5001, 500)) + [np.inf]
LABELS = [f"{i}-{i+499}nt" for i in range(0, 5000, 500)] + ["5000nt+"]

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.size"] = 8
plt.rcParams["svg.fonttype"] = "none"


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=300)
        print(f"Saved: {out}")


print("Loading GTF...")
load_gtf_df(str(GTF_FILE))

print("Loading Salmon quantifications...")
quant_paths = sorted(glob.glob(QUANT_GLOB))
tpm_dfs = []
for path in quant_paths:
    sample_name = Path(path).parent.name.replace("_untrimmed_harmonized", "")
    df = pd.read_csv(path, sep="\t")[["Name", "TPM"]].set_index("Name")
    df.columns = [sample_name]
    tpm_dfs.append(df)

raw_tpm_df = pd.concat(tpm_dfs, axis=1)
lengths_df = pd.read_csv(quant_paths[0], sep="\t")[["Name", "Length"]].set_index("Name")
raw_tpm_df["Length"] = lengths_df["Length"]

print("Filtering to all human transcripts...")
human_ids = [idx for idx in raw_tpm_df.index if "ENST" in idx]
human_tpm = raw_tpm_df.loc[human_ids].copy()
available_cols = [c for c in TOTAL_RNA_COLS if c in human_tpm.columns]
human_tpm_all_filtered = human_tpm[human_tpm[available_cols].mean(axis=1) > 1].copy()
human_tpm_all_filtered["length_category"] = pd.cut(
    human_tpm_all_filtered["Length"],
    bins=BINS,
    labels=LABELS,
)

human_tpm_all_filtered["reference_scaled"] = (
    human_tpm_all_filtered[TOTAL_RNA_REFERENCE].mean(axis=1) *
    (human_tpm_all_filtered["Length"] + 1000)
)
human_tpm_all_filtered["reference_scaled"] = (
    human_tpm_all_filtered["reference_scaled"] * 1e6 /
    human_tpm_all_filtered["reference_scaled"].sum()
)
print(f"  {len(human_tpm_all_filtered)} human transcripts with avg TPM > 1")

fig = plt.figure(figsize=(5, 5))
axes_size_inches = 50 / 25.4
ax = fig.add_axes([0.2, 0.2, axes_size_inches / 5, axes_size_inches / 5])

groups = [
    ("gag:PABP", TOTAL_RNA_GAG_PABP, "#76c7dcff"),
    ("reference", TOTAL_RNA_REFERENCE, "#888888ff"),
    ("gag", TOTAL_RNA_GAG_ONLY, "#ffbc22d8"),
]

for group_name, cols, color in groups:
    available = [c for c in cols if c in human_tpm_all_filtered.columns]
    if not available:
        continue
    means = human_tpm_all_filtered.groupby("length_category", observed=False)[available].mean().mean(axis=1)
    sterr = human_tpm_all_filtered.groupby("length_category", observed=False)[available].sem().mean(axis=1)
    ax.errorbar(
        range(len(LABELS)),
        means,
        yerr=sterr,
        label=group_name,
        color=color,
        marker="o",
        capsize=5,
        markersize=3,
    )
    if group_name == "reference":
        ax.fill_between(range(len(LABELS)), list(means), color=color, alpha=0.2)

ax.set_xticks(range(len(LABELS)))
ax.set_xticklabels(LABELS, rotation=45, ha="right")
ax.set_xlabel("Transcript Length")
ax.set_ylabel("Mean TPM")
ax.legend(frameon=False)
ax.text(0.70, 0.75, f"N={len(human_tpm_all_filtered)}", transform=ax.transAxes, fontsize=8)

save(fig, "suppl5e_total_rna_all_genes_gt1")
plt.close(fig)
print("Done.")
