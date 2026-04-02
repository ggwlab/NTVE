from collections import defaultdict
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from repo_paths import FIGURE2_LOADED_DIR, find_repo_root


def find_root() -> Path:
    return find_repo_root(Path(__file__))


ROOT = find_root()
GENERATED_FIG2_LOADED_DIR = FIGURE2_LOADED_DIR


def get_generated_figure2_loaded_dir() -> Path:
    if not GENERATED_FIG2_LOADED_DIR.exists():
        raise FileNotFoundError(
            f"Missing generated Figure 2 inputs at {GENERATED_FIG2_LOADED_DIR}. "
            "Run `conda run -n ntve python Figure2/figure2_prepare_loaded_data.py` first."
        )
    return GENERATED_FIG2_LOADED_DIR


def create_biotype_dicts(transcript_info_df: pd.DataFrame, gtf_df_data: dict) -> dict:
    gtf_df = gtf_df_data["gtf_df"]
    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype_dict = defaultdict(lambda: "NA", gene_level_gtf.set_index("gene_id")["gene_biotype"].to_dict())
    transcript_to_gene = transcript_info_df.set_index("transcript_id")["gene_id"].to_dict()
    transcript_biotype_dict = {}
    for transcript_id, gene_id in transcript_to_gene.items():
        transcript_biotype_dict[transcript_id] = gene_biotype_dict.get(gene_id, "NA")
    return transcript_biotype_dict


def parse_sample_names_transcript_generic(tpm_df: pd.DataFrame, lims_csv: Path | None = None) -> dict:
    sample_name_mapping = {}
    lims_mapping = {}
    if lims_csv is not None and Path(lims_csv).exists():
        lims_df = pd.read_csv(lims_csv)
        lims_mapping = dict(zip(lims_df["Sample_Name"], lims_df["Sample_NameLIMS"]))

    for col_name in tpm_df.columns:
        col_str = str(col_name)
        if col_str in lims_mapping:
            lims_name = lims_mapping[col_str]
            match = re.match(r"SN_(\d+)_(\d+)", str(lims_name))
            if match:
                sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "SN"}
                continue
            match = re.match(r"L_(\d+)_(\d+)", str(lims_name))
            if match:
                sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "L"}
                continue

        match = re.match(r"SN_S(\d+)_(\d+)", col_str)
        if match:
            sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "SN"}
            continue
        match = re.match(r"Lysate_S(\d+)_(\d+)", col_str)
        if match:
            sample_name_mapping[col_str] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "Lysate"}
    return sample_name_mapping


def organize_samples_by_type_and_number(sample_map: dict, tpm_df: pd.DataFrame) -> dict:
    organized = {}
    for col_name, info in sample_map.items():
        if col_name not in tpm_df.columns:
            continue
        sample_type = info["Type"]
        sample_num = info["Sample"]
        organized.setdefault(sample_type, {}).setdefault(sample_num, []).append(col_name)
    return organized


def create_length_bins() -> tuple[list[float], list[str]]:
    bins = [0, 1000, 2000, 3000, 4000, 5000, float("inf")]
    labels = ["0-999", "1000-1999", "2000-2999", "3000-3999", "4000-4999", "5000+"]
    return bins, labels


def get_sample_pairs(source_samples: list[str], target_samples: list[str], sample_map: dict) -> list[tuple[str, str]]:
    source_by_replicate = defaultdict(list)
    target_by_replicate = defaultdict(list)
    for source_col in source_samples:
        replicate = sample_map.get(source_col, {}).get("Replicate")
        if replicate is not None:
            source_by_replicate[replicate].append(source_col)
    for target_col in target_samples:
        replicate = sample_map.get(target_col, {}).get("Replicate")
        if replicate is not None:
            target_by_replicate[replicate].append(target_col)
    sample_pairs = []
    for replicate in sorted(set(source_by_replicate) & set(target_by_replicate)):
        for source_col in source_by_replicate[replicate]:
            for target_col in target_by_replicate[replicate]:
                sample_pairs.append((source_col, target_col))
    return sample_pairs


def analyze_correlation_by_length(
    tpm_df: pd.DataFrame,
    sample_map: dict,
    transcript_info_df: pd.DataFrame,
    biotype_dict: dict,
    organized: dict,
    source_type: str,
    target_type: str,
    tpm_threshold: float = 1.0,
) -> tuple[dict, list[str]]:
    bins, bin_labels = create_length_bins()
    pc_transcripts = [tid for tid, biotype in biotype_dict.items() if biotype == "protein_coding" and tid in tpm_df.index]
    tpm_pc = tpm_df.loc[pc_transcripts].copy()
    length_dict = transcript_info_df.set_index("transcript_id")["length"].to_dict()
    tpm_pc["length"] = tpm_pc.index.map(lambda tid: length_dict.get(tid, np.nan))
    tpm_pc = tpm_pc.dropna(subset=["length"])
    tpm_pc["length_bin"] = pd.cut(tpm_pc["length"], bins=bins, labels=bin_labels, right=False)

    results = {}
    sample_numbers = sorted(organized.get(source_type, {}))
    for sample_num in sample_numbers:
        source_samples = [s for s in organized[source_type].get(sample_num, []) if s in tpm_pc.columns]
        target_samples = [t for t in organized[target_type].get(sample_num, []) if t in tpm_pc.columns]
        sample_pairs = get_sample_pairs(source_samples, target_samples, sample_map)
        results[sample_num] = {}

        for bin_label in bin_labels:
            bin_data = tpm_pc[tpm_pc["length_bin"] == bin_label].copy()
            if len(bin_data) == 0 or not sample_pairs:
                results[sample_num][bin_label] = {"spearman_r_mean": np.nan, "spearman_r_std": np.nan, "n_transcripts": 0, "spearman_r_values": []}
                continue
            expressed = (bin_data[source_samples] > tpm_threshold).all(axis=1) & (bin_data[target_samples] > tpm_threshold).all(axis=1)
            n_transcripts = int(expressed.sum())
            if n_transcripts <= 1:
                results[sample_num][bin_label] = {"spearman_r_mean": np.nan, "spearman_r_std": np.nan, "n_transcripts": n_transcripts, "spearman_r_values": []}
                continue
            bin_data = bin_data.loc[expressed]
            spearman_values = []
            for source_col, target_col in sample_pairs:
                spearman_r, _ = spearmanr(bin_data[source_col].values, bin_data[target_col].values)
                spearman_values.append(float(spearman_r))
            results[sample_num][bin_label] = {
                "spearman_r_mean": float(np.mean(spearman_values)),
                "spearman_r_std": float(np.std(spearman_values)) if len(spearman_values) > 2 else np.nan,
                "n_transcripts": n_transcripts,
                "spearman_r_values": spearman_values,
            }
    return results, bin_labels


def add_value_dots(ax, x_center: float, values: list[float], width: float) -> None:
    if not values:
        return
    if len(values) == 1:
        x_offsets = np.array([0.0])
    else:
        x_offsets = np.linspace(-width * 0.18, width * 0.18, len(values))
    ax.scatter(x_center + x_offsets, values, s=12, color="black", edgecolors="none", alpha=0.6, zorder=4)
