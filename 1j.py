"""
Figure 1j — Length-resolved transcript abundance (TPM by length bin)
Standalone reproduction of Figure1/length_resolved_analysis_refactored.ipynb
Outputs to: refactoring_roadmap/1j_plots/
"""
import numpy as np
import pandas as pd
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import sys

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "1j_plots"
OUT_DIR.mkdir(exist_ok=True)

GTF_FILE = ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
QUANT_GLOB = str(ROOT / "resources/salmon_harmonized/*/quant.sf.gz")
DOMINANT_ISOFORM = ROOT / "resources/dominant_isoform_genes.csv.gz"

POLYA_REFERENCE_COLS = ['24L006460', '24L006461', '24L006462']
GAG_PABP_SAMPLES    = ['24L006479', '24L006480', '24L006481']
GAG_ONLY_SAMPLES    = ['24L011711', '24L011712', '24L011713']

BINS   = [0] + list(range(500, 5001, 500)) + [np.inf]
LABELS = [f'{i}-{i+499}nt' for i in range(0, 5000, 500)] + ['5000nt+']

# ── Load GTF ──────────────────────────────────────────────────────────────────
print("Loading GTF...")
gtf_dict = load_gtf_df(str(GTF_FILE))
gtf_df = gtf_dict["gtf_df"]
transcript_df = gtf_df.query("feature == 'transcript'").copy().set_index('transcript_id')
transcript_id_to_biotype = transcript_df["transcript_biotype"].to_dict()

# ── Load Salmon quant files ────────────────────────────────────────────────────
quant_paths = sorted(glob.glob(QUANT_GLOB))
print(f"Loading {len(quant_paths)} Salmon quant files...")
tpm_dfs = []
for path in quant_paths:
    sample = path.split("/")[-2].replace("_untrimmed_harmonized", "")
    df = pd.read_csv(path, sep="\t")[['Name', 'TPM']].set_index('Name')
    df.columns = [sample]
    tpm_dfs.append(df)
raw_tpm_df = pd.concat(tpm_dfs, axis=1)
lengths_df = pd.read_csv(quant_paths[0], sep="\t")[['Name', 'Length']].set_index('Name')
raw_tpm_df['Length'] = lengths_df['Length']

# ── Filter: human protein-coding, renormalize ─────────────────────────────────
human_ids = [i for i in raw_tpm_df.index if "ENST" in i]
human_tpm = raw_tpm_df.loc[human_ids].copy()
human_tpm["transcript_biotype"] = [
    transcript_id_to_biotype.get(i.split(".")[0], "unknown") for i in human_tpm.index]
human_tpm_coding = human_tpm.query("transcript_biotype == 'protein_coding'").copy()

sample_cols = [c for c in human_tpm_coding.columns if c != 'Length' and c != 'transcript_biotype']
for col in sample_cols:
    s = human_tpm_coding[col].sum()
    if s > 0:
        human_tpm_coding[col] = human_tpm_coding[col] / s * 1e6

all_cols = list(set(POLYA_REFERENCE_COLS + GAG_PABP_SAMPLES + GAG_ONLY_SAMPLES))
avail = [c for c in all_cols if c in human_tpm_coding.columns]
human_tpm_filtered = human_tpm_coding[human_tpm_coding[avail].mean(axis=1) > 1].copy()
human_tpm_filtered['length_category'] = pd.cut(human_tpm_filtered['Length'], bins=BINS, labels=LABELS)
print(f"  {len(human_tpm_filtered)} protein-coding transcripts (avg TPM>1)")

# ── Plot: polyA SN vs polyA Lysate (Figure 1j) ────────────────────────────────
plt.rcParams.update({'font.size': 8, 'svg.fonttype': 'none'})
axes_size_in = 50 / 25.4

groups = [
    ('gag:PABP', GAG_PABP_SAMPLES, '#76c7dc'),
    ('Reference', POLYA_REFERENCE_COLS, '#888888'),
]

fig = plt.figure(figsize=(5, 5))
ax = fig.add_axes([0.2, 0.2, axes_size_in/5, axes_size_in/5])

for label, cols, color in groups:
    avail_cols = [c for c in cols if c in human_tpm_filtered.columns]
    if not avail_cols:
        continue
    means = human_tpm_filtered.groupby('length_category', observed=False)[avail_cols].mean().mean(axis=1)
    stderr = human_tpm_filtered.groupby('length_category', observed=False)[avail_cols].sem().mean(axis=1)
    ax.errorbar(range(len(LABELS)), means, yerr=stderr,
                label=label, color=color, marker='o', capsize=5, markersize=3)
    if label == 'Reference':
        ax.fill_between(range(len(LABELS)), list(means), color=color, alpha=0.2)

# Pearson r between gag:PABP and reference
pabp_avail = [c for c in GAG_PABP_SAMPLES if c in human_tpm_filtered.columns]
ref_avail  = [c for c in POLYA_REFERENCE_COLS if c in human_tpm_filtered.columns]
pabp_means = human_tpm_filtered.groupby('length_category', observed=False)[pabp_avail].mean().mean(axis=1)
ref_means  = human_tpm_filtered.groupby('length_category', observed=False)[ref_avail].mean().mean(axis=1)
r, p = stats.pearsonr(pabp_means, ref_means)
ax.text(0.05, 0.95, f"Pearson r={r:.3f}\np={p:.2e}", transform=ax.transAxes,
        fontsize=7, va='top')

ax.set_xticks(range(len(LABELS)))
ax.set_xticklabels(LABELS, rotation=45, ha='right', fontsize=7)
ax.set_xlabel('Transcript Length')
ax.set_ylabel('Mean TPM')
ax.legend(frameon=False, fontsize=7)
ax.text(0.70, 0.75, f'N={len(human_tpm_filtered)}', transform=ax.transAxes, fontsize=8)

for ext in ('svg', 'png'):
    out = OUT_DIR / f"1j_length_resolved.{ext}"
    fig.savefig(out, format=ext, bbox_inches='tight', dpi=150)
    print(f"Saved: {out}")
plt.close(fig)
print("Done.")
