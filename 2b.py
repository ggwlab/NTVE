"""
Figure 2b — Benchmarking comparison transcript detection rate curves
Notebook-faithful standalone reproduction of the top row of the combined-line
Elo local-ranking plot from Figure2/rediscovery_analysis_noMT_v3.ipynb.
Outputs to: refactoring_roadmap/2b_plots/
"""

from collections import defaultdict
from pathlib import Path
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent / "2b_plots"
OUT_DIR.mkdir(exist_ok=True)

sys.path.insert(0, str(ROOT))

from ntvetools import load_gtf_df
from refactoring_roadmap_figure2_shared import get_generated_figure2_loaded_dir
from rediscovery_analysis_lib import reorganize_samples_by_number, run_per_sample_analysis

warnings.filterwarnings("ignore")
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 8

name_mapping = {
    "Benchmarking-S1": "NTVE",
    "Benchmarking-S2": "HIV-gag WT",
    "Benchmarking-S3": "TRACE-seq",
    "Benchmarking-S4": "MMLV-gag WT",
    "Benchmarking-S5": "HIV-gag:MCP",
    "Benchmarking-S6": "NLuc",
}


def get_sample_display_name(exp_prefix: str, sample_num: int) -> str:
    key = f"{exp_prefix}-S{sample_num}"
    mapped_name = name_mapping.get(key)
    if mapped_name:
        return f"S{sample_num}: {mapped_name}"
    return f"S{sample_num}"


def aggregate_transcripts_to_genes(tpm_df: pd.DataFrame, transcript_info_df: pd.DataFrame) -> pd.DataFrame:
    transcript_to_gene = transcript_info_df.set_index("transcript_id")["gene_id"].to_dict()
    tpm_with_gene = tpm_df.copy()
    tpm_with_gene["gene_id"] = tpm_with_gene.index.map(lambda tid: transcript_to_gene.get(tid, ""))
    return tpm_with_gene.groupby("gene_id").sum(numeric_only=True)


def create_gene_df(gene_tpm: pd.DataFrame, gene_numreads: pd.DataFrame, gtf_df_data: dict, gene_id_to_gene_name: dict, gene_id_is_mt: dict) -> pd.DataFrame:
    gtf_df = gtf_df_data["gtf_df"]
    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype_dict = defaultdict(lambda: "NA", gene_level_gtf.set_index("gene_id")["gene_biotype"].to_dict())
    gene_df = gene_tpm.copy()
    gene_df.index.name = "GeneID"
    gene_df["GeneName"] = gene_df.index.map(lambda gid: gene_id_to_gene_name.get(gid, gid))
    gene_df["GeneBiotype"] = gene_df.index.map(lambda gid: gene_biotype_dict.get(gid, "NA"))
    gene_df["is_MT"] = gene_df.index.map(lambda gid: gene_id_is_mt.get(gid, False))
    for col in gene_numreads.columns:
        gene_df[col.replace("_NumReads", "") + "_NumReads"] = gene_numreads[col]
    return gene_df


def filter_genes(gene_df: pd.DataFrame) -> pd.DataFrame:
    custom_genes = {"mGreenLantern", "Gag_common_core", "HIV_gag_MCP", "MMLV_gag", "minigag_PABP", "TetOn3G-V9I_mScarlet", "eUnaG2"}
    df = gene_df.copy()
    df = df[df["GeneBiotype"] == "protein_coding"].copy()
    df = df[df["is_MT"] == False].copy()
    df = df[~df["GeneName"].isin(custom_genes)].copy()
    return df


def renormalize_tpm(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    numeric_cols = [col for col in out.columns if col not in ["GeneName", "GeneBiotype", "is_MT"] and not col.endswith("_NumReads")]
    for col in numeric_cols:
        total = out[col].sum()
        if total > 0:
            out[col] = out[col] * 1e6 / total
    return out


def parse_sample_names(numeric_cols_elo: list[str]) -> dict:
    sample_name_mapping = {}
    for col in numeric_cols_elo:
        match = re.match(r"SN_S(\d+)_(\d+)", str(col))
        if match:
            sample_name_mapping[col] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "SN", "Experiment": "Benchmarking comparison", "ColumnName": col}
            continue
        match = re.match(r"Lysate_S(\d+)_(\d+)", str(col))
        if match:
            sample_name_mapping[col] = {"Sample": int(match.group(1)), "Replicate": int(match.group(2)), "Type": "Lysate", "Experiment": "Benchmarking comparison", "ColumnName": col}
    return sample_name_mapping


def save_plot(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
        print(f"Saved: {out}")


def plot_top_row_local_ranking(results, analysis_types, colors, exp_prefix, linestyles, linewidths):
    fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)

    for analysis_idx, (source_type, target_type, direction) in enumerate(analysis_types[:2]):
        ax = axes[analysis_idx]
        analysis_key = f"{source_type}→{target_type}_{direction}"

        for sample_num in sorted(results.keys()):
            result = results[sample_num][analysis_key]
            line_label = f"{get_sample_display_name(exp_prefix, sample_num)} (AUC={result['auc']:.2f})"
            ax.plot(
                result["quantiles"],
                result["rates"],
                label=line_label,
                color=colors[sample_num],
                linewidth=linewidths[sample_num],
                linestyle=linestyles[sample_num],
                alpha=0.95,
            )

        ax.axhline(y=50, color="gray", linestyle="--", alpha=0.4, linewidth=1)
        ax.set_xlim(-5, 105)
        ax.set_ylim(-5, 105)
        ax.grid(True, alpha=0.2)
        ax.set_title(f"{direction.upper()}: {source_type}→{target_type}", fontweight="bold", fontsize=8)
        ax.set_xlabel("Expression Quantile (%)", fontweight="bold", fontsize=8)
        if analysis_idx == 0:
            ax.set_ylabel("Rediscovery Rate (%)", fontweight="bold", fontsize=8)
        ax.legend(loc="lower right", fontsize=8, frameon=True, framealpha=0.85)

    fig.suptitle(
        "Benchmarking comparison - Local Ranking Rediscovery",
        fontsize=8,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()
    return fig


print("Loading GTF and transcript matrices...")
gtf_df_data = load_gtf_df(str(ROOT / "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"))
gene_id_to_gene_name = defaultdict(str, gtf_df_data["gene_id_to_name"])
gene_id_is_mt = gtf_df_data.get("gene_id_is_mt", {})
output_dir = get_generated_figure2_loaded_dir()
benchmark_tpm = pd.read_csv(output_dir / "elo_transcript_tpm_matrix.csv", index_col=0)
benchmark_numreads = pd.read_csv(output_dir / "elo_transcript_numreads_matrix.csv", index_col=0)
benchmark_transcript_info = pd.read_csv(output_dir / "elo_transcript_info.csv")
benchmark_gene_tpm = aggregate_transcripts_to_genes(benchmark_tpm, benchmark_transcript_info)
benchmark_gene_numreads = aggregate_transcripts_to_genes(benchmark_numreads, benchmark_transcript_info)
benchmark_gene_df = create_gene_df(benchmark_gene_tpm, benchmark_gene_numreads, gtf_df_data, gene_id_to_gene_name, gene_id_is_mt)
gene_df_elo_pc = renormalize_tpm(filter_genes(benchmark_gene_df))

numeric_cols = [col for col in gene_df_elo_pc.columns if col not in ["GeneName", "GeneBiotype", "is_MT"] and not col.endswith("_NumReads")]
sample_name_mapping = parse_sample_names(numeric_cols)
elo_sn_samples = [col for col, meta in sample_name_mapping.items() if meta["Type"] == "SN"]
elo_lysate_samples = [col for col, meta in sample_name_mapping.items() if meta["Type"] == "Lysate"]
elo_organized = {
    "SN": reorganize_samples_by_number(elo_sn_samples, sample_name_mapping),
    "Lysate": reorganize_samples_by_number(elo_lysate_samples, sample_name_mapping),
}

elo_colors = {
    1: "#77C8DC",
    2: "#BDC7E3",
    3: "#ECB0A3",
    4: "#D7E6D5",
    5: "#B4C7B0",
    6: "black",
}
elo_linestyles = {
    1: "-",
    2: (0, (8, 3)),
    3: (0, (3, 2)),
    4: "-.",
    5: (0, (8, 3, 1, 3, 1, 3)),
    6: (0, (8, 3)),
}
elo_linewidths = {
    1: 3.5,
    2: 2.0,
    3: 2.0,
    4: 2.0,
    5: 2.0,
    6: 2.0,
}
analysis_types_elo = [("Lysate", "SN", "bottom"), ("Lysate", "SN", "top"), ("SN", "Lysate", "bottom"), ("SN", "Lysate", "top")]

print("Running benchmarking comparison per-sample rediscovery analysis (local ranking)...")
elo_local_results = run_per_sample_analysis(elo_organized, gene_df_elo_pc, sample_name_mapping, "Lysate", "SN")
fig = plot_top_row_local_ranking(elo_local_results, analysis_types_elo, elo_colors, "Benchmarking", elo_linestyles, elo_linewidths)
save_plot(fig, "2b_benchmarking_local_ranking_toprow")
plt.close(fig)
print("Done.")
