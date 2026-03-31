"""
Figure 7 supporting 2vs2 causal DESeq2 + fgsea pipeline.

Recreates the intermediate and downstream data from:
- Figure7/DESeq2_Timepoint_Comparison_GSEA_Window-2vs2.ipynb

Outputs under Figure7/deseq2_2vs2_causal_output/:
- dayX/deseq2_results.csv
- dayX/deseq2_results_significant.csv
- dayX/normalized_counts.csv
- gene_rankings_padj1.0.csv
- fgsea_analysis/output/fgsea_c2_all_timepoints.csv
- fgsea_analysis/output/fgsea_c5_all_timepoints.csv
- fgsea_analysis/output/summary_per_timepoint.csv
- fgsea_analysis/output/top_c5_pathways_per_timepoint.csv
- fgsea_analysis/output/pathways_per_timepoint.pdf
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from fig7_shared import ROOT, build_cardio_count_matrix

WINDOW_SIZE = 4
PADJ_THRESHOLD = 1.0
MIN_GENE_SET_SIZE = 15
MAX_GENE_SET_SIZE = 500
TIMEPOINTS = list(range(1, 10))

TEMP_BASE = ROOT / "Figure7" / "deseq2_temp"
OUTPUT_DIR = ROOT / "Figure7" / "deseq2_2vs2_causal_output"
FGSEA_INPUT_DIR = OUTPUT_DIR / "fgsea_analysis" / "input"
FGSEA_OUTPUT_DIR = OUTPUT_DIR / "fgsea_analysis" / "output"
DESEQ_DOCKER_DIR = ROOT / "resources" / "docker_DeSeq_with_Shrinkage"
FGSEA_DOCKER_DIR = ROOT / "resources" / "docker_fgsea"
DESEQ_IMAGE = "deseq2-with-shrinkage:latest"
FGSEA_IMAGE = "fgsea:latest"


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def create_deseq2_r_script_timepoint_comparison(timepoint: int, before_days: list[int]) -> str:
    before_label = f"{before_days[0]}-{before_days[-1]}" if len(before_days) > 1 else str(before_days[0])
    return f"""#!/usr/bin/env Rscript
library(DESeq2)

cat("\\n=== DESeq2 Analysis: Day {timepoint} vs Days {before_label} ===\\n")
cat("Method: Design ~ Condition with apeglm shrinkage\\n\\n")

setwd("/workspace/input")
output_dir <- "/workspace/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count_matrix <- read.csv("count_matrix.csv", row.names = 1, check.names = FALSE)
sample_metadata <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)
sample_metadata <- sample_metadata[order(match(sample_metadata$Sample, colnames(count_matrix))), ]
count_matrix <- count_matrix[, sample_metadata$Sample]

rownames(sample_metadata) <- sample_metadata$Sample
sample_metadata$Condition <- factor(sample_metadata$Condition, levels = c("Before", "Current"))

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData = sample_metadata,
  design = ~ Condition
)
dds <- DESeq(dds)

coef_name <- resultsNames(dds)[2]
res <- results(dds, name = coef_name, alpha = 0.05)

tryCatch({{
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
}}, error = function(e) {{
  res_shrunk <<- res
}})

gene_id_to_name_csv <- read.csv("gene_id_to_name.csv", row.names = 1)
gene_id_to_name <- setNames(gene_id_to_name_csv$GeneName, rownames(gene_id_to_name_csv))

results_df <- as.data.frame(res_shrunk)
results_df$Gene <- rownames(results_df)
results_df$GeneName <- gene_id_to_name[results_df$Gene]
results_df <- results_df[, c("Gene", "GeneName", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

write.csv(results_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)
sig_genes <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
write.csv(sig_genes, file.path(output_dir, "deseq2_results_significant.csv"), row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- data.frame(Gene = rownames(norm_counts), norm_counts)
write.csv(norm_counts_df, file.path(output_dir, "normalized_counts.csv"), row.names = FALSE)
"""


def create_fgsea_r_script() -> str:
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

run_fgsea_timepoint <- function(rnk_file, pathways, pathway_name, min_size, max_size) {{
  ranks_df <- fread(rnk_file)
  ranks <- setNames(ranks_df$rank, ranks_df$gene)
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size)
  tp <- sub(".*day([0-9]+)\\\\.rnk$", "\\\\1", rnk_file)
  res$timepoint <- as.integer(tp)
  res$collection <- pathway_name
  res
}}

rnk_files <- list.files("/workspace/input", pattern = "^ranked_genes_day[0-9]+\\\\.rnk$", full.names = TRUE)
rnk_files <- sort(rnk_files)

c2_res <- rbindlist(lapply(rnk_files, run_fgsea_timepoint, pathways = c2_pathways, pathway_name = "C2", min_size = MIN_SIZE, max_size = MAX_SIZE), fill = TRUE)
c5_res <- rbindlist(lapply(rnk_files, run_fgsea_timepoint, pathways = c5_pathways, pathway_name = "C5", min_size = MIN_SIZE, max_size = MAX_SIZE), fill = TRUE)

fwrite(c2_res, file.path(output_dir, "fgsea_c2_all_timepoints.csv"), sep = "\\t")
fwrite(c5_res, file.path(output_dir, "fgsea_c5_all_timepoints.csv"), sep = "\\t")
"""


def determine_window(tp: int) -> tuple[list[int], list[int]]:
    if tp == 1:
        return [1], [0]
    if tp == 2:
        return [2], [0, 1]
    return [tp - 1, tp], [tp - 3, tp - 2]


def prepare_per_timepoint_inputs(count_pc: pd.DataFrame, full_metadata: pd.DataFrame, gene_id_to_name: dict) -> dict[int, pd.DataFrame]:
    all_results: dict[int, pd.DataFrame] = {}
    for tp in TIMEPOINTS:
        current_days, before_days = determine_window(tp)
        rows = []
        for _, row in full_metadata.iterrows():
            if row["Day"] in before_days:
                condition = "Before"
            elif row["Day"] in current_days:
                condition = "Current"
            else:
                continue
            rows.append(
                {
                    "Sample": row["Sample"],
                    "FullName": row["FullName"],
                    "Day": row["Day"],
                    "Replicate": row["Replicate"],
                    "Condition": condition,
                }
            )

        sample_metadata = pd.DataFrame(rows).sort_values(["Condition", "Day", "Replicate"]).reset_index(drop=True)
        count_matrix = count_pc[sample_metadata["FullName"].values].copy()
        count_matrix.columns = sample_metadata["Sample"].values

        temp_dir = TEMP_BASE / f"day{tp}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        count_matrix.to_csv(temp_dir / "count_matrix.csv")
        sample_metadata[["Sample", "Condition"]].to_csv(temp_dir / "sample_metadata.csv", index=False)
        pd.DataFrame(
            {"GeneName": [gene_id_to_name.get(gene_id, gene_id) for gene_id in count_matrix.index]},
            index=count_matrix.index,
        ).to_csv(temp_dir / "gene_id_to_name.csv")
        (temp_dir / "run_deseq2.R").write_text(create_deseq2_r_script_timepoint_comparison(tp, before_days))
    return all_results


def rank_results(results_by_timepoint: dict[int, pd.DataFrame]) -> pd.DataFrame:
    all_data = []
    for tp in TIMEPOINTS:
        if tp not in results_by_timepoint:
            continue
        results_df = results_by_timepoint[tp].dropna(subset=["log2FoldChange", "padj"]).copy()
        filtered = results_df.copy() if PADJ_THRESHOLD >= 1.0 else results_df[results_df["padj"] < PADJ_THRESHOLD].copy()
        upregulated = filtered[filtered["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
        downregulated = filtered[filtered["log2FoldChange"] < 0].sort_values("log2FoldChange", ascending=True)
        for direction, subset in [("Upregulated", upregulated), ("Downregulated", downregulated)]:
            if subset.empty:
                continue
            out = subset[["Gene", "GeneName", "log2FoldChange", "padj"]].copy()
            out["Rank"] = range(1, len(out) + 1)
            out["Direction"] = direction
            out["Timepoint"] = tp
            all_data.append(out[["Timepoint", "Direction", "Rank", "Gene", "GeneName", "log2FoldChange", "padj"]])
    combined = pd.concat(all_data, ignore_index=True)
    combined.to_csv(OUTPUT_DIR / f"gene_rankings_padj{PADJ_THRESHOLD}.csv", index=False)
    return combined


def export_fgsea_rank_files(rankings: pd.DataFrame) -> None:
    FGSEA_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    FGSEA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for name in ["c2.all.v2025.1.Hs.json", "c5.go.v2025.1.Hs.json"]:
        candidates = [
            ROOT / "resources" / name,
            ROOT / "Figure7" / "deseq2_2vs2_causal_output" / "fgsea_analysis" / "input" / name,
            ROOT / "Suppl21" / "fgsea_analysis" / "input" / name,
        ]
        src = next((path for path in candidates if path.exists()), None)
        if src is None:
            raise FileNotFoundError(f"Missing gene-set file: tried {candidates}")
        dst = FGSEA_INPUT_DIR / name
        if src.resolve() != dst.resolve():
            shutil.copy(src, dst)

    for tp in TIMEPOINTS:
        tp_df = rankings[rankings["Timepoint"] == tp].copy()
        if tp_df.empty:
            continue
        tp_df["rank_value"] = tp_df["log2FoldChange"]
        ranked = tp_df[["GeneName", "rank_value"]].rename(columns={"GeneName": "gene", "rank_value": "rank"})
        ranked["abs_rank"] = ranked["rank"].abs()
        ranked = ranked.sort_values("abs_rank", ascending=False).drop_duplicates(subset="gene", keep="first")
        ranked = ranked[["gene", "rank"]].sort_values("rank", ascending=False)
        ranked.to_csv(FGSEA_INPUT_DIR / f"ranked_genes_day{tp}.rnk", sep="\t", index=False)

    (FGSEA_INPUT_DIR / "run_fgsea.R").write_text(create_fgsea_r_script())


def summarise_fgsea() -> None:
    c2 = pd.read_csv(FGSEA_OUTPUT_DIR / "fgsea_c2_all_timepoints.csv", sep="\t")
    c5 = pd.read_csv(FGSEA_OUTPUT_DIR / "fgsea_c5_all_timepoints.csv", sep="\t")

    summary_rows = []
    top_rows = []
    for tp in sorted(c5["timepoint"].dropna().unique()):
        tp_int = int(tp)
        c2_tp = c2[c2["timepoint"] == tp_int]
        c5_tp = c5[c5["timepoint"] == tp_int]
        summary_rows.append(
            {
                "Timepoint": f"Day {tp_int}",
                "C2_Significant": int((c2_tp["padj"] < 0.05).sum()),
                "C2_Total": int(len(c2_tp)),
                "C5_Significant": int((c5_tp["padj"] < 0.05).sum()),
                "C5_Total": int(len(c5_tp)),
            }
        )
        c5_sig = c5_tp[c5_tp["padj"] < 0.05].sort_values("pval")
        if not c5_sig.empty:
            top = c5_sig.iloc[0]
            top_rows.append(
                {
                    "Timepoint": f"Day {tp_int}",
                    "pathway": top["pathway"],
                    "pval": top["pval"],
                    "padj": top["padj"],
                    "ES": top["ES"],
                    "NES": top["NES"],
                    "size": top["size"],
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(FGSEA_OUTPUT_DIR / "summary_per_timepoint.csv", index=False)
    pd.DataFrame(top_rows).to_csv(FGSEA_OUTPUT_DIR / "top_c5_pathways_per_timepoint.csv", index=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    c2_sig_counts = c2[c2["padj"] < 0.05].groupby("timepoint").size().reset_index(name="count")
    c5_sig_counts = c5[c5["padj"] < 0.05].groupby("timepoint").size().reset_index(name="count")
    ax1.bar(c2_sig_counts["timepoint"], c2_sig_counts["count"], color="steelblue", edgecolor="black")
    ax1.set_xlabel("Timepoint (Days)", fontweight="bold")
    ax1.set_ylabel("Significant Pathways (padj < 0.05)", fontweight="bold")
    ax1.set_title(f"C2 Pathways Enriched Per Timepoint\n(Window Size = {WINDOW_SIZE})", fontweight="bold")
    ax1.grid(axis="y", alpha=0.3)
    ax2.bar(c5_sig_counts["timepoint"], c5_sig_counts["count"], color="coral", edgecolor="black")
    ax2.set_xlabel("Timepoint (Days)", fontweight="bold")
    ax2.set_ylabel("Significant Pathways (padj < 0.05)", fontweight="bold")
    ax2.set_title(f"C5 (GO) Pathways Enriched Per Timepoint\n(Window Size = {WINDOW_SIZE})", fontweight="bold")
    ax2.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FGSEA_OUTPUT_DIR / "pathways_per_timepoint.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TEMP_BASE.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Preparing cardio count matrix from quantfiles...")
    count_matrix, full_metadata, ref = build_cardio_count_matrix(include_day10=False)
    count_pc_full = pd.read_csv(ROOT / "Figure7" / "cardio_deseq2_cubicsplines_input" / "count_matrix.csv", index_col=0) if (ROOT / "Figure7" / "cardio_deseq2_cubicsplines_input" / "count_matrix.csv").exists() else None
    full_name_order = full_metadata["FullName"].values
    count_pc = pd.DataFrame(index=count_matrix.index)
    for short_name, full_name in zip(count_matrix.columns, full_name_order):
        count_pc[full_name] = count_matrix[short_name].values

    prepare_per_timepoint_inputs(count_pc, full_metadata, ref["gene_id_to_name"])

    print("Building DESeq2 and fgsea Docker images...")
    run(["docker", "build", "-t", DESEQ_IMAGE, "."], cwd=DESEQ_DOCKER_DIR)
    run(["docker", "build", "-t", FGSEA_IMAGE, "."], cwd=FGSEA_DOCKER_DIR)

    results_by_timepoint: dict[int, pd.DataFrame] = {}
    for tp in TIMEPOINTS:
        temp_dir = TEMP_BASE / f"day{tp}"
        tp_output_dir = OUTPUT_DIR / f"day{tp}"
        tp_output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Running DESeq2 for day {tp}...")
        run(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{temp_dir}:/workspace/input:rw",
                "-v",
                f"{tp_output_dir}:/workspace/output:rw",
                DESEQ_IMAGE,
                "Rscript",
                "/workspace/input/run_deseq2.R",
            ]
        )
        results_by_timepoint[tp] = pd.read_csv(tp_output_dir / "deseq2_results.csv")

    rankings = rank_results(results_by_timepoint)
    export_fgsea_rank_files(rankings)

    print("Running fgsea...")
    run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{FGSEA_INPUT_DIR}:/workspace/input:rw",
            "-v",
            f"{FGSEA_OUTPUT_DIR}:/workspace/output:rw",
            FGSEA_IMAGE,
            "Rscript",
            "/workspace/input/run_fgsea.R",
        ]
    )
    summarise_fgsea()
    print("Done.")


if __name__ == "__main__":
    main()
