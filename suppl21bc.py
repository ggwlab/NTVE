"""
Supplementary Figure 21b,c — fgsea per timepoint.

Refactors:
- Suppl21/GSEA_fgsea_per_timepoint.ipynb

Default behaviour:
1. Ensure rankings exist (via suppl21bc_rankings.py)
2. Export .rnk files
3. Run fgsea in Docker
4. Recreate the notebook's intermediate output tables
5. Emit the notebook's supporting significance-count barplots

Primary notebook-style outputs:
- refactoring_roadmap/suppl21bc_output/fgsea_c2_all_timepoints.csv
- refactoring_roadmap/suppl21bc_output/fgsea_c5_all_timepoints.csv
- refactoring_roadmap/suppl21bc_output/c2_aggregated_across_timepoints.csv
- refactoring_roadmap/suppl21bc_output/c5_aggregated_across_timepoints.csv
- refactoring_roadmap/suppl21bc_output/c2_fisher_combined_pvalues.csv
- refactoring_roadmap/suppl21bc_output/c5_fisher_combined_pvalues.csv
- refactoring_roadmap/suppl21bc_output/c2_stouffer_combined_pvalues.csv
- refactoring_roadmap/suppl21bc_output/c5_stouffer_combined_pvalues.csv

Manuscript-facing plots:
- refactoring_roadmap/suppl21bc_plots/suppl21b_go_day_resolved.{png,svg}
- refactoring_roadmap/suppl21bc_plots/suppl21c_go_whole_trajectory.{png,svg}
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from ntvetools import load_gtf_df
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests

def find_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here.parent, here.parent.parent):
        if (candidate / "resources").exists():
            return candidate
    return here.parent.parent


def first_existing_path(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


ROOT = find_root()
WORK_DIR = Path(__file__).parent
OUT_DATA = WORK_DIR / "suppl21bc_output"
OUT_PLOTS = WORK_DIR / "suppl21bc_plots"
FGSEA_INPUT = OUT_DATA / "fgsea_input"
FGSEA_OUTPUT = OUT_DATA / "fgsea_output"
OUT_DATA.mkdir(exist_ok=True)
OUT_PLOTS.mkdir(exist_ok=True)
FGSEA_INPUT.mkdir(exist_ok=True)
FGSEA_OUTPUT.mkdir(exist_ok=True)

PADJ_THRESHOLD = 0.5
MIN_GENE_SET_SIZE = 15
MAX_GENE_SET_SIZE = 500
FGSEA_IMAGE = "fgsea:latest"
FGSEA_DOCKER_DIR = ROOT / "resources" / "docker_fgsea"
TOP_N_PER_DAY = 3

RANKINGS_CSV = OUT_DATA / "all_gene_rankings_log2_centered.csv"
DESEQ2_CSV = ROOT / "Figure7" / "cardio_deseq2_cubicsplines_output" / "deseq2_results.csv"
GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

PANEL_B_ROWS = [
    ("contractile muscle fiber", "GOCC_CONTRACTILE_MUSCLE_FIBER"),
    ("I-band", "GOCC_I_BAND"),
    ("muscle system process", "GOBP_MUSCLE_SYSTEM_PROCESS"),
    ("muscle contraction", "GOBP_MUSCLE_CONTRACTION"),
    ("actin filament based movement", "GOBP_ACTIN_FILAMENT_BASED_MOVEMENT"),
    ("cardiac muscle contractions", "GOBP_CARDIAC_MUSCLE_CONTRACTION"),
    ("supramolecular complex", "GOCC_SUPRAMOLECULAR_COMPLEX"),
    ("supramolecular polymer", "GOCC_SUPRAMOLECULAR_POLYMER"),
    ("cellular assembly in morphogenesis", "GOBP_CELLULAR_COMPONENT_ASSEMBLY_INVOLVED_IN_MORPHOGENESIS"),
]

PANEL_C_LABELS = {
    "GOBP_MUSCLE_STRUCTURE_DEVELOPMENT": "muscle structure development",
    "GOBP_MUSCLE_ORGAN_DEVELOPMENT": "muscle organ development",
    "GOBP_CARDIAC_MUSCLE_TISSUE_DEVELOPMENT": "cardiac muscle tissue development",
    "GOBP_MUSCLE_TISSUE_DEVELOPMENT": "muscle tissue development",
    "GOBP_CARDIAC_MUSCLE_CONTRACTION": "cardiac muscle contractions",
    "GOBP_STRIATED_MUSCLE_CONTRACTION": "striated muscle contractions",
    "GOBP_HEART_PROCESS": "heart process",
    "GOBP_ACTIN_FILAMENT_BASED_MOVEMENT": "actin filament based movement",
    "GOBP_MUSCLE_SYSTEM_PROCESS": "muscle system process",
    "GOBP_HEART_MORPHOGENESIS": "heart morphogenesis",
    "GOBP_CARDIAC_VENTRICLE_DEVELOPMENT": "cardiac ventricle development",
    "GOCC_SUPRAMOLECULAR_POLYMER": "supramolecular polymer",
    "GOBP_CIRCULATORY_SYSTEM_PROCESS": "circulatory system process",
    "GOCC_I_BAND": "I-band",
    "GOCC_CONTRACTILE_MUSCLE_FIBER": "contractile muscle fiber",
}


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "svg"):
        out = OUT_PLOTS / f"{stem}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


def ensure_rankings() -> None:
    if not RANKINGS_CSV.exists():
        run([sys.executable, str(WORK_DIR / "suppl21bc_rankings.py")])


def export_rnk_files() -> None:
    print("Loading rankings and DESeq2 results...")
    if not DESEQ2_CSV.exists():
        raise FileNotFoundError(
            "Missing cardio spline DESeq2 results. Run `fig7c_pipeline.py` first to generate "
            f"{ROOT / 'Figure7' / 'cardio_deseq2_cubicsplines_output' / 'deseq2_results.csv'}."
        )
    rankings = pd.read_csv(RANKINGS_CSV)
    deseq2_results = pd.read_csv(DESEQ2_CSV)

    if PADJ_THRESHOLD < 1.0:
        significant_genes = deseq2_results[deseq2_results["padj"] < PADJ_THRESHOLD]["Gene"].unique()
        rankings = rankings[rankings["Gene"].isin(significant_genes)].copy()
        print(f"Filtered rankings to {rankings['Gene'].nunique()} genes with padj < {PADJ_THRESHOLD}")

    gtf_data = load_gtf_df(str(GTF_PATH))
    gene_id_to_symbol = gtf_data["gene_id_to_name"]
    rankings["GeneSymbol"] = rankings["Gene"].map(gene_id_to_symbol)
    rankings = rankings.dropna(subset=["GeneSymbol"]).copy()

    timepoints = sorted(rankings["Timepoint"].unique())
    print(f"Exporting ranked lists for {len(timepoints)} timepoints...")
    for tp in timepoints:
        tp_data = rankings[rankings["Timepoint"] == tp].copy()
        tp_up = tp_data[tp_data["Direction"] == "Upregulated"].sort_values("Log2-Centered", ascending=False)
        tp_down = tp_data[tp_data["Direction"] == "Downregulated"].sort_values("Log2-Centered", ascending=True)
        tp_combined = pd.concat([tp_up, tp_down])
        ranked_list = tp_combined[["GeneSymbol", "Log2-Centered"]].copy()
        ranked_list.columns = ["gene", "rank"]
        ranked_list["abs_rank"] = ranked_list["rank"].abs()
        ranked_list = ranked_list.sort_values("abs_rank", ascending=False).drop_duplicates(subset="gene", keep="first")
        ranked_list = ranked_list[["gene", "rank"]].sort_values("rank", ascending=False)
        out = FGSEA_INPUT / f"ranked_genes_day{tp}.rnk"
        ranked_list.to_csv(out, sep="\t", index=False)
        print(f"  Day {tp}: {len(ranked_list)} genes → {out.name}")

    for name in ("c2.all.v2025.1.Hs.json", "c5.go.v2025.1.Hs.json"):
        src = ROOT / "resources" / name
        if not src.exists():
            raise FileNotFoundError(f"Missing gene-set file: {src}")
        shutil.copy(src, FGSEA_INPUT / name)


def build_r_script() -> str:
    return f"""#!/usr/bin/env Rscript
suppressPackageStartupMessages({{
  library(fgsea)
  library(data.table)
  library(BiocParallel)
  library(jsonlite)
}})

register(MulticoreParam(workers = max(1, parallel::detectCores() - 1)))

setwd("/workspace/input")
output_dir <- "/workspace/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_SIZE <- {MIN_GENE_SET_SIZE}
MAX_SIZE <- {MAX_GENE_SET_SIZE}

c2_json <- fromJSON("/workspace/input/c2.all.v2025.1.Hs.json")
c2_pathways <- lapply(names(c2_json), function(name) c2_json[[name]]$geneSymbols)
names(c2_pathways) <- names(c2_json)

c5_json <- fromJSON("/workspace/input/c5.go.v2025.1.Hs.json")
c5_pathways <- lapply(names(c5_json), function(name) c5_json[[name]]$geneSymbols)
names(c5_pathways) <- names(c5_json)

run_fgsea_timepoint <- function(rnk_file, pathways) {{
  timepoint <- gsub("ranked_genes_day(\\\\d+)\\\\.rnk", "\\\\1", basename(rnk_file))
  ranks_df <- read.table(rnk_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
  ranks <- setNames(ranks_df$rank, ranks_df$gene)
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = ranks, minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0.0)
  res <- as.data.frame(res)
  res <- res[order(res$pval), ]
  res$timepoint <- as.integer(timepoint)
  res
}}

rnk_files <- sort(list.files(pattern="^ranked_genes_day.*\\\\.rnk$", full.names=TRUE))
c2_results_list <- lapply(rnk_files, function(f) run_fgsea_timepoint(f, c2_pathways))
c5_results_list <- lapply(rnk_files, function(f) run_fgsea_timepoint(f, c5_pathways))
c2_all <- do.call(rbind, c2_results_list)
c5_all <- do.call(rbind, c5_results_list)
fwrite(c2_all, file.path(output_dir, "fgsea_c2_all_timepoints.csv"), sep="\\t")
fwrite(c5_all, file.path(output_dir, "fgsea_c5_all_timepoints.csv"), sep="\\t")
"""


def run_fgsea() -> None:
    run(["docker", "build", "-t", FGSEA_IMAGE, "."], cwd=FGSEA_DOCKER_DIR)
    r_script_path = FGSEA_INPUT / "run_fgsea.R"
    r_script_path.write_text(build_r_script())
    run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{FGSEA_INPUT}:/workspace/input:rw",
            "-v",
            f"{FGSEA_OUTPUT}:/workspace/output:rw",
            FGSEA_IMAGE,
            "Rscript",
            "/workspace/input/run_fgsea.R",
        ]
    )


def plot_outputs() -> None:
    c5 = pd.read_csv(FGSEA_OUTPUT / "fgsea_c5_all_timepoints.csv", sep="\t")
    c5.to_csv(OUT_DATA / "fgsea_c5_all_timepoints.csv", index=False)
    pd.read_csv(FGSEA_OUTPUT / "fgsea_c2_all_timepoints.csv", sep="\t").to_csv(
        OUT_DATA / "fgsea_c2_all_timepoints.csv", index=False
    )
    plot_panel_b(c5)
    plot_panel_c()


def plot_panel_b(c5: pd.DataFrame) -> None:
    row_labels = [label for label, _ in PANEL_B_ROWS]
    pathway_to_label = {pathway: label for label, pathway in PANEL_B_ROWS}
    days = list(range(10))

    fig, ax = plt.subplots(figsize=(8.4, 3.4))
    ax.set_xlim(-0.6, len(days) - 0.4)
    ax.set_ylim(-0.5, len(row_labels) - 0.5)
    ax.invert_yaxis()

    for y in range(len(row_labels)):
        ax.plot([-0.45, len(days) - 0.55], [y, y], color="#efefef", lw=14, solid_capstyle="round", zorder=0)
    for x in range(len(days) + 1):
        ax.axvline(x - 0.5, color="#999999", lw=0.6, ls=(0, (2, 2)), zorder=0)

    selected_rows = []
    for day in days:
        day_hits = (
            c5[(c5["timepoint"] == day) & (c5["padj"] < 0.05)]
            .sort_values("pval")
            .head(TOP_N_PER_DAY)
        )
        for rank_within_day, (_, row) in enumerate(day_hits.iterrows(), start=1):
            pathway = row["pathway"]
            if pathway not in pathway_to_label:
                continue
            x = int(row["timepoint"])
            y = row_labels.index(pathway_to_label[pathway])
            symbol = "+" if row["NES"] > 0 else "−"
            ax.text(
                x,
                y,
                symbol,
                ha="center",
                va="center",
                fontsize=12,
                bbox=dict(boxstyle="circle,pad=0.12", fc="white", ec="black", lw=1),
                zorder=3,
            )
            selected_rows.append(
                {
                    "day": day,
                    "rank_within_day": rank_within_day,
                    "pathway": pathway,
                    "label": pathway_to_label[pathway],
                    "pval": row["pval"],
                    "padj": row["padj"],
                    "NES": row["NES"],
                }
            )

    ax.set_xticks(range(len(days)))
    ax.set_xticklabels([f"d{d}" for d in days], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.tick_params(length=0)
    ax.set_title("Timepoint (days)", fontsize=10, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.tight_layout()
    save(fig, "suppl21b_go_day_resolved")
    plt.close(fig)

    selected_df = pd.DataFrame(selected_rows)
    selected_out = OUT_DATA / "suppl21b_selected_hits_top3.csv"
    selected_df.to_csv(selected_out, index=False)
    print(f"Saved: {selected_out}")


def plot_panel_c() -> None:
    agg = pd.read_csv(OUT_DATA / "c5_aggregated_across_timepoints.csv")
    top15 = agg.sort_values("sum_rank").head(15).copy()
    top15["label"] = top15["pathway"].map(PANEL_C_LABELS).fillna(
        top15["pathway"]
        .str.replace("GOBP_", "", regex=False)
        .str.replace("GOCC_", "", regex=False)
        .str.replace("_", " ", regex=False)
        .str.lower()
    )
    top15 = top15.iloc[::-1]

    fig, ax = plt.subplots(figsize=(4.2, 5.2))
    ax.barh(top15["label"], top15["sum_rank"], color="#cfcfcf", edgecolor="none", height=0.72)
    ax.scatter(top15["sum_rank"], top15["label"], color="black", s=8, zorder=3)
    ax.set_xlim(0, 820)
    ax.set_xticks([0, 400, 800])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("rank sum", fontsize=10, fontweight="bold", labelpad=6)
    ax.set_title("Gene ontology enrichment analysis\nof whole trajectory (n=10480 pathways)", fontsize=10, pad=28)
    ax.tick_params(axis="y", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1)

    fig.tight_layout()
    save(fig, "suppl21c_go_whole_trajectory")
    plt.close(fig)


def build_aggregated_table(df: pd.DataFrame) -> pd.DataFrame:
    clean = df.dropna(subset=["pval"]).copy()
    clean["rank_in_timepoint"] = clean.groupby("timepoint")["pval"].rank(method="min")
    aggregated = (
        clean.groupby("pathway")
        .agg(
            sum_pval=("pval", "sum"),
            sum_rank=("rank_in_timepoint", "sum"),
            n_timepoints=("timepoint", "count"),
            n_significant=("padj", lambda x: (x < 0.05).sum()),
            mean_NES=("NES", "mean"),
        )
        .reset_index()
        .sort_values(["sum_rank", "sum_pval"])
    )
    aggregated["rank_by_sum_pval"] = aggregated["sum_pval"].rank(method="min")
    aggregated["rank_by_sum_rank"] = aggregated["sum_rank"].rank(method="min")
    return aggregated[
        [
            "pathway",
            "sum_pval",
            "sum_rank",
            "n_timepoints",
            "n_significant",
            "mean_NES",
            "rank_by_sum_pval",
            "rank_by_sum_rank",
        ]
    ]


def build_combined_pvalue_table(df: pd.DataFrame, method: str) -> pd.DataFrame:
    clean = df.dropna(subset=["pval"]).copy()
    rows = []
    for pathway, group in clean.groupby("pathway"):
        pvals = group["pval"].to_numpy()
        fisher_pval = combine_pvalues(pvals, method="fisher")[1]
        stouffer_pval = combine_pvalues(pvals, method="stouffer")[1]
        rows.append(
            {
                "pathway": pathway,
                "fisher_pval": fisher_pval,
                "stouffer_pval": stouffer_pval,
                "min_pval": pvals.min(),
                "n_timepoints": len(pvals),
                "n_significant": int((group["padj"] < 0.05).sum()),
                "mean_NES": group["NES"].mean(),
            }
        )

    combined = pd.DataFrame(rows)
    padj_col = "fisher_padj" if method == "fisher" else "stouffer_padj"
    source_col = "fisher_pval" if method == "fisher" else "stouffer_pval"
    combined[padj_col] = multipletests(combined[source_col], method="fdr_bh")[1]
    combined = combined.sort_values(source_col)

    if method == "fisher":
        return combined[
            [
                "pathway",
                "fisher_pval",
                "stouffer_pval",
                "min_pval",
                "n_timepoints",
                "n_significant",
                "mean_NES",
                "fisher_padj",
            ]
        ]

    return combined[
        [
            "pathway",
            "stouffer_pval",
            "min_pval",
            "n_timepoints",
            "n_significant",
            "mean_NES",
            "stouffer_padj",
        ]
    ]


def export_notebook_intermediates() -> None:
    print("Computing notebook-style aggregated and combined-pvalue tables...")
    c2 = pd.read_csv(FGSEA_OUTPUT / "fgsea_c2_all_timepoints.csv", sep="\t")
    c5 = pd.read_csv(FGSEA_OUTPUT / "fgsea_c5_all_timepoints.csv", sep="\t")

    c2.to_csv(OUT_DATA / "fgsea_c2_all_timepoints.csv", index=False)
    c5.to_csv(OUT_DATA / "fgsea_c5_all_timepoints.csv", index=False)

    outputs = {
        "c2_aggregated_across_timepoints.csv": build_aggregated_table(c2),
        "c5_aggregated_across_timepoints.csv": build_aggregated_table(c5),
        "c2_fisher_combined_pvalues.csv": build_combined_pvalue_table(c2, method="fisher"),
        "c5_fisher_combined_pvalues.csv": build_combined_pvalue_table(c5, method="fisher"),
        "c2_stouffer_combined_pvalues.csv": build_combined_pvalue_table(c2, method="stouffer"),
        "c5_stouffer_combined_pvalues.csv": build_combined_pvalue_table(c5, method="stouffer"),
    }

    for name, df in outputs.items():
        out = OUT_DATA / name
        df.to_csv(out, index=False)
        print(f"Saved: {out}")


def main() -> None:
    ensure_rankings()
    export_rnk_files()
    run_fgsea()
    export_notebook_intermediates()
    plot_outputs()
    print("Done.")


if __name__ == "__main__":
    main()
