"""
Generate the Figure 2 / Supplementary Figure 11 cached quantification matrices
directly from primary quantification inputs.

This replaces the old dependency on notebook-exported:
  - Figure2/loaded_data
  - Suppl11/loaded_data

Outputs to:
  refactoring_roadmap/figure2_loaded_data/
"""

from collections import defaultdict
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df

OUT_DIR = Path(__file__).parent / "figure2_loaded_data"
OUT_DIR.mkdir(exist_ok=True)

QUANT_DIR = ROOT / "resources" / "quantfiles_filtered_pipeline" / "homo_sapiens"
LIMS_CSV = ROOT / "resources" / "Project_1716_lims_simplified_Spam2_without_missing.csv"
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


def load_mappings() -> tuple[defaultdict, defaultdict, defaultdict, defaultdict]:
    print("Loading GTF reference data...")
    gtf_data = load_gtf_df(str(GTF_CSV))
    gtf_df = gtf_data["gtf_df"]

    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype = gene_level_gtf.set_index("gene_id")["gene_biotype"].to_dict()
    transcript_to_gene = gtf_df.set_index("transcript_id")["gene_id"].to_dict()
    transcript_to_gene_name = gtf_df.set_index("transcript_id")["gene_name"].to_dict()
    gene_to_name = gtf_data["gene_id_to_name"]

    for gene_id in SYNTHETIC_GENES:
        transcript_to_gene[gene_id] = gene_id
        transcript_to_gene_name[gene_id] = gene_id
        gene_biotype[gene_id] = "protein_coding"
        gene_to_name[gene_id] = gene_id

    print(f"Loaded mappings for {len(transcript_to_gene)} transcripts")
    return (
        defaultdict(str, transcript_to_gene),
        defaultdict(str, transcript_to_gene_name),
        defaultdict(str, gene_to_name),
        defaultdict(lambda: "NA", gene_biotype),
    )


def load_quant_df(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def aggregate_sample_to_gene(sample_file: Path, transcript_to_gene: dict) -> tuple[dict, dict]:
    df = load_quant_df(sample_file).copy()
    df["GeneID"] = df["Name"].str.split(".").str[0].map(lambda tid: transcript_to_gene[tid])
    agg = df.groupby("GeneID")[["TPM", "NumReads"]].sum()
    return agg["TPM"].to_dict(), agg["NumReads"].to_dict()


def load_transcript_sample(sample_file: Path) -> pd.DataFrame:
    df = load_quant_df(sample_file).copy()
    df["transcript_id"] = df["Name"].str.split(".").str[0]
    return df[["transcript_id", "Length", "TPM", "NumReads"]].rename(
        columns={"Length": "length", "TPM": "tpm", "NumReads": "numreads"}
    )


def build_project1716_matrices(
    transcript_to_gene: dict,
    gene_to_name: dict,
    gene_biotype: dict,
) -> None:
    print("Building Project 1716 matrices...")
    sample_metadata_csv = pd.read_csv(LIMS_CSV)
    sample_names = set(sample_metadata_csv["Sample_Name"].unique())
    sample_files = sorted(
        path for path in (QUANT_DIR / f"{name}.sf.gz" for name in sample_names) if path.exists()
    )
    print(f"Found {len(sample_files)} Project 1716 quantification files")

    sample_tpm_data: dict[str, dict] = {}
    sample_numreads_data: dict[str, dict] = {}
    sample_metadata: dict[str, dict] = {}
    all_genes: set[str] = set()

    for sf_path in sample_files:
        sample_name = sf_path.name.replace(".sf.gz", "")
        csv_row = sample_metadata_csv[sample_metadata_csv["Sample_Name"] == sample_name]
        if csv_row.empty:
            continue
        tpm_dict, numreads_dict = aggregate_sample_to_gene(sf_path, transcript_to_gene)
        sample_tpm_data[sample_name] = tpm_dict
        sample_numreads_data[sample_name] = numreads_dict
        sample_metadata[sample_name] = csv_row.iloc[0].to_dict()
        all_genes.update(g for g in tpm_dict if g)

    all_genes = sorted(all_genes)
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
    tpm_matrix_pc["GeneName"] = [gene_to_name.get(gene_id, gene_id) for gene_id in tpm_matrix_pc.index]
    numreads_pc = numreads_matrix[is_protein_coding].copy()

    all_transcript_dfs = []
    for sf_path in sample_files:
        sample_name = sf_path.name.replace(".sf.gz", "")
        df = load_quant_df(sf_path).copy()
        df["transcript_id"] = df["Name"].str.split(".").str[0]
        df["sample_name"] = sample_name
        all_transcript_dfs.append(df[["transcript_id", "Length", "TPM", "NumReads", "sample_name"]])

    all_transcripts_combined = pd.concat(all_transcript_dfs, ignore_index=True)
    transcript_info_df = all_transcripts_combined[["transcript_id", "Length"]].drop_duplicates().reset_index(drop=True)
    transcript_info_df.columns = ["transcript_id", "length"]
    transcript_info_df["gene_id"] = transcript_info_df["transcript_id"].map(lambda tid: transcript_to_gene.get(tid, ""))
    transcript_info_df["gene_name"] = transcript_info_df["gene_id"].map(
        lambda gid: "" if gid in SYNTHETIC_GENES else gene_to_name.get(gid, "")
    )

    transcript_tpm_matrix = all_transcripts_combined.pivot_table(
        index="transcript_id", columns="sample_name", values="TPM", fill_value=0
    )
    transcript_numreads_matrix = all_transcripts_combined.pivot_table(
        index="transcript_id", columns="sample_name", values="NumReads", fill_value=0
    )

    sample_meta_df = pd.DataFrame.from_dict(sample_metadata, orient="index")
    sample_meta_df.index.name = "sample_name"
    gene_info = pd.DataFrame(
        {
            "GeneID": all_genes,
            "GeneName": [gene_to_name.get(g, g) for g in all_genes],
            "GeneBiotype": [gene_biotype.get(g, "NA") for g in all_genes],
        }
    )

    tpm_matrix_pc.to_csv(OUT_DIR / "tpm_matrix_protein_coding.csv")
    numreads_pc.to_csv(OUT_DIR / "numreads_matrix_protein_coding.csv")
    sample_meta_df.to_csv(OUT_DIR / "sample_metadata.csv")
    gene_info.to_csv(OUT_DIR / "gene_info.csv", index=False)
    transcript_info_df.to_csv(OUT_DIR / "transcript_info.csv", index=False)
    transcript_tpm_matrix.to_csv(OUT_DIR / "transcript_tpm_matrix.csv")
    transcript_numreads_matrix.to_csv(OUT_DIR / "transcript_numreads_matrix.csv")
    print(f"Saved Project 1716 matrices to {OUT_DIR}")


def parse_elo_sample_name(name: str) -> dict | None:
    if name.startswith("SN_S"):
        sample_type = "SN"
        suffix = name[4:]
    elif name.startswith("Lysate_S"):
        sample_type = "Lysate"
        suffix = name[8:]
    else:
        return None

    parts = suffix.split("_")
    try:
        construct = int(parts[0])
        replicate = int(parts[1])
    except (IndexError, ValueError):
        return None

    if construct < 1 or construct > 6 or replicate < 1 or replicate > 3:
        return None
    return {"sample_type": sample_type, "construct": construct, "replicate": replicate}


def build_elo_matrices(
    transcript_to_gene: dict,
    gene_to_name: dict,
    gene_biotype: dict,
) -> None:
    print("Building Elo benchmarking matrices...")
    elo_files = sorted(QUANT_DIR.glob("SN_S*.sf.gz")) + sorted(QUANT_DIR.glob("Lysate_S*.sf.gz"))

    elo_tpm_data: dict[str, dict] = {}
    elo_numreads_data: dict[str, dict] = {}
    elo_metadata: dict[str, dict] = {}
    elo_genes: set[str] = set()

    for sf_path in elo_files:
        sample_name = sf_path.name.replace(".sf.gz", "")
        meta = parse_elo_sample_name(sample_name)
        if meta is None:
            continue
        tpm_dict, numreads_dict = aggregate_sample_to_gene(sf_path, transcript_to_gene)
        elo_tpm_data[sample_name] = tpm_dict
        elo_numreads_data[sample_name] = numreads_dict
        elo_metadata[sample_name] = meta
        elo_genes.update(g for g in tpm_dict if g)

    elo_genes = sorted(elo_genes)
    elo_tpm_matrix = pd.DataFrame(
        {
            sample_name: [elo_tpm_data[sample_name].get(gene_id, 0) for gene_id in elo_genes]
            for sample_name in elo_tpm_data
        },
        index=elo_genes,
    )
    elo_numreads_matrix = pd.DataFrame(
        {
            sample_name: [elo_numreads_data[sample_name].get(gene_id, 0) for gene_id in elo_genes]
            for sample_name in elo_numreads_data
        },
        index=elo_genes,
    )

    is_protein_coding = [gene_biotype.get(gene_id, "NA") == "protein_coding" for gene_id in elo_tpm_matrix.index]
    elo_tpm_matrix_pc = elo_tpm_matrix[is_protein_coding].copy()
    elo_tpm_matrix_pc["GeneName"] = [gene_to_name.get(gene_id, gene_id) for gene_id in elo_tpm_matrix_pc.index]
    elo_numreads_pc = elo_numreads_matrix[is_protein_coding].copy()

    elo_all_transcript_dfs = []
    for sf_path in elo_files:
        sample_name = sf_path.name.replace(".sf.gz", "")
        if parse_elo_sample_name(sample_name) is None:
            continue
        df = load_quant_df(sf_path).copy()
        df["transcript_id"] = df["Name"].str.split(".").str[0]
        df["sample_name"] = sample_name
        elo_all_transcript_dfs.append(df[["transcript_id", "Length", "TPM", "NumReads", "sample_name"]])

    elo_all_transcripts_combined = pd.concat(elo_all_transcript_dfs, ignore_index=True)
    elo_transcript_info_df = elo_all_transcripts_combined[["transcript_id", "Length"]].drop_duplicates().reset_index(drop=True)
    elo_transcript_info_df.columns = ["transcript_id", "length"]
    elo_transcript_info_df["gene_id"] = elo_transcript_info_df["transcript_id"].map(lambda tid: transcript_to_gene.get(tid, ""))
    elo_transcript_info_df["gene_name"] = elo_transcript_info_df["gene_id"].map(
        lambda gid: "" if gid in SYNTHETIC_GENES else gene_to_name.get(gid, "")
    )

    elo_transcript_tpm_matrix = elo_all_transcripts_combined.pivot_table(
        index="transcript_id", columns="sample_name", values="TPM", fill_value=0
    )
    elo_transcript_numreads_matrix = elo_all_transcripts_combined.pivot_table(
        index="transcript_id", columns="sample_name", values="NumReads", fill_value=0
    )

    elo_meta_df = pd.DataFrame.from_dict(elo_metadata, orient="index")
    elo_meta_df.index.name = "sample_name"

    elo_tpm_matrix_pc.to_csv(OUT_DIR / "elo_tpm_matrix_protein_coding.csv")
    elo_numreads_pc.to_csv(OUT_DIR / "elo_numreads_matrix_protein_coding.csv")
    elo_meta_df.to_csv(OUT_DIR / "elo_sample_metadata.csv")
    elo_transcript_info_df.to_csv(OUT_DIR / "elo_transcript_info.csv", index=False)
    elo_transcript_tpm_matrix.to_csv(OUT_DIR / "elo_transcript_tpm_matrix.csv")
    elo_transcript_numreads_matrix.to_csv(OUT_DIR / "elo_transcript_numreads_matrix.csv")
    print(f"Saved Elo matrices to {OUT_DIR}")


def main() -> None:
    transcript_to_gene, transcript_to_gene_name, gene_to_name, gene_biotype = load_mappings()
    build_project1716_matrices(transcript_to_gene, gene_to_name, gene_biotype)
    build_elo_matrices(transcript_to_gene, gene_to_name, gene_biotype)
    print("Done.")


if __name__ == "__main__":
    main()
