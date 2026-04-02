"""
Rebuild the Figure 6 trilineage pickle directly from raw Salmon quantification
files, including the SN_IPSC baseline samples that the original notebook
discovered but failed to parse into the saved metadata.

Outputs:
  refactoring_roadmap/figure6_three_lineage_data.pkl

Saved object structure intentionally mirrors Figure6/pipeline/three_lineage_data.pkl
so downstream scripts can consume it with minimal changes.
"""

from collections import defaultdict
from pathlib import Path
import pickle
import sys

import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df

OUT_PKL = Path(__file__).parent / "figure6_three_lineage_data.pkl"
QUANTIF_DIR = ROOT / "resources" / "quantfiles_filtered_pipeline" / "homo_sapiens"
GTF_CSV = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

SYNTHETIC_GENES = [
    "TetOn3G-V9I_mScarlet",
    "mGreenLantern",
    "Gag_common_core",
    "MMLV_gag",
    "HIV_gag_MCP",
    "minigag_PABP",
    "eUnaG2",
]


def load_mappings() -> tuple[defaultdict, defaultdict, defaultdict, dict]:
    print("Loading GTF annotation data...")
    gtf_data = load_gtf_df(str(GTF_CSV))
    gtf_df = gtf_data["gtf_df"]

    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype = gene_level_gtf.set_index("gene_id")["gene_biotype"].to_dict()
    transcript_to_gene = gtf_df.set_index("transcript_id")["gene_id"].to_dict()
    transcript_to_gene_name = gtf_df.set_index("transcript_id")["gene_name"].to_dict()
    gene_id_to_name = gtf_data["gene_id_to_name"]

    for gene_id in SYNTHETIC_GENES:
        transcript_to_gene[gene_id] = gene_id
        transcript_to_gene_name[gene_id] = gene_id
        gene_biotype[gene_id] = "protein_coding"
        gene_id_to_name[gene_id] = gene_id

    print(f"Loaded {len(transcript_to_gene)} transcript mappings")
    return (
        defaultdict(str, transcript_to_gene),
        defaultdict(str, transcript_to_gene_name),
        defaultdict(str, gene_id_to_name),
        gene_biotype,
    )


def parse_sample_name(filename: str) -> dict | None:
    """
    Parse both trilineage differentiation SN samples and undifferentiated
    trilineage baseline iPSC supernatant samples.

    Accepted forms:
      SN_ec_d1_6.sf.gz     -> lineage='ectoderm', day=1, clone=6
      SN_IPSC_6.sf.gz      -> lineage='ipsc', day=0, clone=6
    """
    name = filename.replace(".sf.gz", "").replace(".sf", "")

    if not name.startswith("SN_"):
        return None
    if "_dox" in name:
        return None

    parts = name.split("_")

    if len(parts) == 3 and parts[1].upper() == "IPSC":
        try:
            clone = int(parts[2])
        except ValueError:
            return None
        return {"lineage": "ipsc", "day": 0, "clone": clone}

    if len(parts) < 4:
        return None

    lineage_short = parts[1]
    lineage_map = {"mes": "mesoderm", "ec": "ectoderm", "en": "endoderm"}
    if lineage_short not in lineage_map:
        return None

    day_str = parts[2]
    if not day_str.startswith("d"):
        return None

    try:
        day = int(day_str[1:])
        clone = int(parts[3])
    except ValueError:
        return None

    return {"lineage": lineage_map[lineage_short], "day": day, "clone": clone}


def load_quant_df(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_and_aggregate_sample(
    sample_file_path: Path,
    transcript_to_gene_map: dict,
) -> tuple[dict, dict]:
    df = load_quant_df(sample_file_path).copy()

    tpm_dict: dict[str, float] = {}
    numreads_dict: dict[str, float] = {}

    for gene_name in SYNTHETIC_GENES:
        if gene_name in df["Name"].values:
            row = df[df["Name"] == gene_name].iloc[0]
            tpm_dict[gene_name] = float(row["TPM"])
            numreads_dict[gene_name] = float(row["NumReads"])

    transcript_ids = df["Name"].str.split(".").str[0]
    df["GeneID"] = transcript_ids.map(lambda tid: transcript_to_gene_map[tid])
    agg = df.groupby("GeneID")[["TPM", "NumReads"]].sum()

    tpm_dict.update(agg["TPM"].to_dict())
    numreads_dict.update(agg["NumReads"].to_dict())
    return tpm_dict, numreads_dict


def main() -> None:
    transcript_to_gene, transcript_to_gene_name, gene_id_to_name, gene_biotype = load_mappings()

    print(f"Quantification directory: {QUANTIF_DIR}")
    sample_files = sorted(QUANTIF_DIR.glob("SN_*.sf.gz"))
    print(f"Found {len(sample_files)} potential SN samples")

    sample_tpm_data: dict[str, dict] = {}
    sample_numreads_data: dict[str, dict] = {}
    sample_metadata: dict[str, dict] = {}
    all_genes: set[str] = set()

    skipped = 0
    loaded = 0
    for sf_path in sample_files:
        sample_name = sf_path.name.replace(".sf.gz", "")
        meta = parse_sample_name(sample_name)
        if meta is None:
            skipped += 1
            continue

        try:
            tpm_dict, numreads_dict = load_and_aggregate_sample(sf_path, transcript_to_gene)
        except Exception as exc:
            print(f"  Warning: failed to load {sample_name}: {exc}")
            skipped += 1
            continue

        sample_tpm_data[sample_name] = tpm_dict
        sample_numreads_data[sample_name] = numreads_dict
        sample_metadata[sample_name] = meta
        all_genes.update(g for g in tpm_dict if g)
        loaded += 1

    print(f"Loaded {loaded} samples, skipped {skipped}")
    all_genes = sorted(all_genes)
    print(f"Total unique genes: {len(all_genes)}")

    tpm_matrix = pd.DataFrame(
        {
            sample_name: [sample_tpm_data[sample_name].get(gene_id, 0) for gene_id in all_genes]
            for sample_name in sample_tpm_data
        },
        index=all_genes,
    )
    numreads_matrix = pd.DataFrame(
        {
            sample_name: [sample_numreads_data[sample_name].get(gene_id, 0) for gene_id in all_genes]
            for sample_name in sample_numreads_data
        },
        index=all_genes,
    )

    is_protein_coding = [gene_biotype.get(gene_id, "NA") == "protein_coding" for gene_id in tpm_matrix.index]
    tpm_matrix_pc = tpm_matrix[is_protein_coding].copy()
    numreads_matrix_pc = numreads_matrix[is_protein_coding].copy()

    for gene_id in SYNTHETIC_GENES:
        if gene_id not in tpm_matrix_pc.index:
            tpm_matrix_pc.loc[gene_id] = 0.0
            numreads_matrix_pc.loc[gene_id] = 0.0

    cpm_matrix_pc = (numreads_matrix_pc / numreads_matrix_pc.sum(axis=0)) * 1e6

    tpm_matrix_pc["GeneName"] = [gene_id_to_name.get(gene_id, gene_id) for gene_id in tpm_matrix_pc.index]
    numreads_matrix_pc["GeneName"] = [gene_id_to_name.get(gene_id, gene_id) for gene_id in numreads_matrix_pc.index]
    cpm_matrix_pc["GeneName"] = [gene_id_to_name.get(gene_id, gene_id) for gene_id in cpm_matrix_pc.index]

    sample_meta_df = pd.DataFrame.from_dict(sample_metadata, orient="index").sort_values(
        by=["lineage", "clone", "day"]
    )
    sample_meta_df.index.name = "sample_name"

    print(f"Sample metadata shape: {sample_meta_df.shape}")
    print(f"Lineages: {sorted(sample_meta_df['lineage'].unique().tolist())}")
    print(f"Days: {sorted(sample_meta_df['day'].unique().tolist())}")
    print(f"Clones: {sorted(sample_meta_df['clone'].unique().tolist())}")

    data_to_save = {
        "tpm_matrix_pc": tpm_matrix_pc,
        "sample_meta_df": sample_meta_df,
        "gene_biotype_dict": gene_biotype,
        "transcript_id_to_gene_id": dict(transcript_to_gene),
        "transcript_id_to_gene_name": dict(transcript_to_gene_name),
        "gene_id_to_gene_name": dict(gene_id_to_name),
        "gtf_mappings": {
            "transcript_id_to_gene_id": dict(transcript_to_gene),
            "transcript_id_to_gene_name": dict(transcript_to_gene_name),
            "gene_biotype": gene_biotype,
        },
        "numreads_matrix_pc": numreads_matrix_pc,
        "cpm_matrix_pc": cpm_matrix_pc,
    }

    with OUT_PKL.open("wb") as f:
        pickle.dump(data_to_save, f)

    print(f"Saved: {OUT_PKL}")
    print(f"  - tpm_matrix_pc: {tpm_matrix_pc.shape}")
    print(f"  - numreads_matrix_pc: {numreads_matrix_pc.shape}")
    print(f"  - cpm_matrix_pc: {cpm_matrix_pc.shape}")
    print(f"  - sample_meta_df: {sample_meta_df.shape}")


if __name__ == "__main__":
    main()
