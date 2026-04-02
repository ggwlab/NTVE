"""
Supplementary Figure 5d — Total-RNA coding TPM by transcript length
Caption-driven standalone reproduction of the total-RNA coding plot block from
Figure1/length_resolved_analysis_refactored.ipynb.

Outputs to: refactoring_roadmap/suppl5d_plots/
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
sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "suppl5d_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "suppl5d_csv"
CSV_DIR.mkdir(exist_ok=True)

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


def renormalize_tpm(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = df.copy()
    for col in columns:
        if col in result.columns:
            total = result[col].sum()
            if total > 0:
                result[col] = (result[col] / total) * 1e6
    return result


print("Loading GTF...")
gtf_df = load_gtf_df(str(GTF_FILE))["gtf_df"]
transcript_df = gtf_df.query("feature == 'transcript'").copy().set_index("transcript_id")
transcript_id_to_biotype = transcript_df["transcript_biotype"].to_dict()

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

print("Filtering to human protein-coding transcripts...")
human_ids = [idx for idx in raw_tpm_df.index if "ENST" in idx]
human_tpm = raw_tpm_df.loc[human_ids].copy()
human_tpm["transcript_biotype"] = [
    transcript_id_to_biotype.get(idx.split(".")[0], "unknown")
    for idx in human_tpm.index
]
human_tpm_coding = human_tpm.query("transcript_biotype == 'protein_coding'").copy()
human_tpm_coding = renormalize_tpm(human_tpm_coding, list(human_tpm_coding.columns[:-2]))

available_cols = [c for c in TOTAL_RNA_COLS if c in human_tpm_coding.columns]
human_tpm_filtered = human_tpm_coding[human_tpm_coding[available_cols].mean(axis=1) > 1].copy()
human_tpm_filtered["length_category"] = pd.cut(human_tpm_filtered["Length"], bins=BINS, labels=LABELS)
print(f"  {len(human_tpm_filtered)} coding transcripts with avg TPM > 1")

fig = plt.figure(figsize=(5, 5))
axes_size_inches = 50 / 25.4
ax = fig.add_axes([0.2, 0.2, axes_size_inches / 5, axes_size_inches / 5])

groups = [
    ("gag:PABP", TOTAL_RNA_GAG_PABP, "#76c7dcff"),
    ("reference", TOTAL_RNA_REFERENCE, "#888888ff"),
    ("gag", TOTAL_RNA_GAG_ONLY, "#ffbc22d8"),
]
summary_rows = []

for group_name, cols, color in groups:
    available = [c for c in cols if c in human_tpm_filtered.columns]
    if not available:
        continue
    means = human_tpm_filtered.groupby("length_category", observed=False)[available].mean().mean(axis=1)
    stds = human_tpm_filtered.groupby("length_category", observed=False)[available].mean().std(axis=1, ddof=1)
    sterr = human_tpm_filtered.groupby("length_category", observed=False)[available].sem().mean(axis=1)
    summary_rows.extend(
        {
            "group": group_name,
            "length_category": category,
            "n_replicates": len(available),
            "mean_tpm": mean,
            "std_tpm": std,
            "stderr_tpm": err,
        }
        for category, mean, std, err in zip(LABELS, means, stds, sterr)
    )
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
ax.text(0.70, 0.75, f"N={len(human_tpm_filtered)}", transform=ax.transAxes, fontsize=8)
pd.DataFrame(summary_rows).to_csv(CSV_DIR / "suppl5d_length_resolved_summary.csv", index=False)

save(fig, "suppl5d_total_rna_coding_gt1")
plt.close(fig)
print("Done.")
