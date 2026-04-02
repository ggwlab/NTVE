from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ntvetools import load_gtf_df

def find_root() -> Path:
    return Path(__file__).resolve().parent.parent


ROOT = find_root()
QUANT_DIR = ROOT / "resources" / "quantfiles_filtered_pipeline" / "homo_sapiens"
GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

CUSTOM_GENES = {
    "TetOn3G-V9I_mScarlet": "TetOn3G-V9I_mScarlet",
    "mGreenLantern": "mGreenLantern",
    "Gag_common_core": "Gag_common_core",
    "MMLV_gag": "MMLV_gag",
    "HIV_gag_MCP": "HIV_gag_MCP",
    "minigag_PABP": "minigag_PABP",
    "eUnaG2": "eUnaG2",
}

_CM_RE = re.compile(r"^cm_(\d+)_(\d+)")


def parse_cm_name(name: str) -> tuple[int, int] | None:
    match = _CM_RE.match(name)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def load_cardio_reference() -> dict:
    gtf_data = load_gtf_df(str(GTF_PATH))
    gtf_df = gtf_data["gtf_df"]

    transcript_to_gene_id = (
        gtf_df[gtf_df["feature"] == "transcript"]
        .dropna(subset=["transcript_id"])
        .set_index("transcript_id")["gene_id"]
        .to_dict()
    )
    transcript_to_gene_name = (
        gtf_df[gtf_df["feature"] == "transcript"]
        .dropna(subset=["transcript_id"])
        .set_index("transcript_id")["gene_name"]
        .to_dict()
    )
    gene_biotype_dict = (
        gtf_df[gtf_df["feature"] == "gene"][["gene_id", "gene_biotype"]]
        .drop_duplicates()
        .set_index("gene_id")["gene_biotype"]
        .to_dict()
    )

    for custom_id, custom_name in CUSTOM_GENES.items():
        transcript_to_gene_id[custom_id] = custom_id
        transcript_to_gene_name[custom_id] = custom_name
        gene_biotype_dict[custom_id] = "custom"

    return {
        "gtf_data": gtf_data,
        "transcript_to_gene_id": defaultdict(str, transcript_to_gene_id),
        "transcript_to_gene_name": defaultdict(str, transcript_to_gene_name),
        "gene_biotype": defaultdict(lambda: "NA", gene_biotype_dict),
        "gene_id_to_name": gtf_data["gene_id_to_name"],
    }


def load_and_process_sample(
    sample_path: Path,
    transcript_to_gene_id: dict,
) -> pd.Series:
    df = pd.read_csv(sample_path, sep="\t")
    df["GeneID"] = [transcript_to_gene_id[t.split(".")[0]] for t in df["Name"]]
    gene_level = df.groupby("GeneID")[["NumReads"]].sum()
    return gene_level["NumReads"]


def build_cardio_count_matrix(include_day10: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ref = load_cardio_reference()
    sample_files = sorted(QUANT_DIR.glob("cm_*.sf.gz"))
    print(f"Found {len(sample_files)} cm_* quantfiles")

    all_genes: set[str] = set()
    sample_data: dict[str, pd.Series] = {}
    sample_meta: dict[str, tuple[int, int]] = {}

    for sf in sample_files:
        sample_name = sf.stem.replace(".sf", "")
        meta = parse_cm_name(sample_name)
        if meta is None:
            continue
        day, _ = meta
        if day == 10 and not include_day10:
            continue
        counts = load_and_process_sample(sf, ref["transcript_to_gene_id"])
        sample_data[sample_name] = counts
        sample_meta[sample_name] = meta
        all_genes.update(counts.index.tolist())

    all_genes_sorted = sorted(g for g in all_genes if g)
    count_df = pd.DataFrame(
        {sample: sample_data[sample].reindex(all_genes_sorted, fill_value=0) for sample in sample_data},
        index=all_genes_sorted,
    ).astype(int)

    keep_mask = [ref["gene_biotype"][gene_id] in {"protein_coding", "custom"} for gene_id in count_df.index]
    count_pc = count_df.loc[keep_mask].copy()

    metadata_rows = []
    for full_name in count_pc.columns:
        day, rep = sample_meta[full_name]
        metadata_rows.append(
            {
                "Sample": f"cm_{day}_{rep}",
                "FullName": full_name,
                "Day": day,
                "Time": day,
                "Replicate": rep,
            }
        )

    sample_metadata = pd.DataFrame(metadata_rows).sort_values(["Day", "Replicate"]).reset_index(drop=True)
    count_matrix = count_pc[sample_metadata["FullName"].values].copy()
    count_matrix.columns = sample_metadata["Sample"].values
    count_matrix.index.name = "Gene"

    return count_matrix, sample_metadata, ref
