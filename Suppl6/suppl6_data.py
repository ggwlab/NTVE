"""
Suppl6 data generation — aggregates cm_* Salmon quantfiles to a
protein-coding, non-MT gene-level TPM matrix and saves it to
suppl6_plots/suppl6_cm_tpm_matrix.csv.

Run this once before suppl6ab.py.
"""

import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df

QUANT_DIR = ROOT / "resources" / "quantfiles_filtered_pipeline" / "homo_sapiens"
GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
OUT_DIR = Path(__file__).parent / "suppl6_plots"
OUT_DIR.mkdir(exist_ok=True)
OUT_FILE = OUT_DIR / "suppl6_cm_tpm_matrix.csv"

# Sequencing-run suffix pattern: _MKDL... or similar all-caps 4-letter code
_SUFFIX_RE = re.compile(r"_[A-Z]{4}\d+.*\.sf\.gz$")


def short_name(path: Path) -> str:
    """Strip sequencing-run suffix to recover the biological sample name."""
    return _SUFFIX_RE.sub("", path.name)


# --- Reference data ---
print("Loading GTF data...")
gtf_data = load_gtf_df(str(GTF_PATH))
gtf_df = gtf_data["gtf_df"]

gene_meta = gtf_df.query("feature == 'gene'")[["gene_id", "gene_biotype", "seqname"]]
protein_coding_ids = set(gene_meta.query("gene_biotype == 'protein_coding'")["gene_id"])
mt_gene_ids = set(gene_meta.query("seqname == 'MT'")["gene_id"])
keep_genes = protein_coding_ids - mt_gene_ids
print(f"Protein-coding non-MT gene IDs: {len(keep_genes)}")

# transcript → gene mapping (strip version suffix)
transcript_to_gene = (
    gtf_df.dropna(subset=["transcript_id"])
    .set_index("transcript_id")["gene_id"]
    .to_dict()
)

gene_id_to_name = gtf_data["gene_id_to_name"]

# --- Load quantfiles ---
cm_files = sorted(QUANT_DIR.glob("cm_*.sf.gz"))
print(f"Found {len(cm_files)} cm_* quantfiles")

tpm_by_sample: dict[str, pd.Series] = {}
for sf in cm_files:
    name = short_name(sf)
    df = pd.read_csv(sf, sep="\t")
    # Strip transcript version (ENST00000xxx.2 → ENST00000xxx)
    df["gene_id"] = df["Name"].str.split(".").str[0].map(transcript_to_gene)
    df = df.dropna(subset=["gene_id"])
    gene_tpm = df.groupby("gene_id")["TPM"].sum()
    tpm_by_sample[name] = gene_tpm
    print(f"  {name}: {len(gene_tpm)} genes")

# --- Assemble matrix, filter, annotate ---
tpm_matrix = pd.DataFrame(tpm_by_sample).fillna(0.0)
tpm_matrix = tpm_matrix.loc[tpm_matrix.index.isin(keep_genes)].copy()
tpm_matrix.insert(0, "GeneName", [gene_id_to_name.get(g, g) for g in tpm_matrix.index])
tpm_matrix.index.name = "gene_id"

tpm_matrix.to_csv(OUT_FILE)
print(f"\nSaved: {OUT_FILE}")
print(f"Matrix shape: {tpm_matrix.shape[0]} genes × {tpm_matrix.shape[1] - 1} samples")
