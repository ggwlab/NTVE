"""
Figure 4c-e upstream DESeq2 generation.

Recreates the per-sample TP2 vs TP3 DESeq2 Wald results from:
  Figure4/neuron_TP2_vs_TP3_DESeq2-Wald.ipynb

Outputs to:
  refactoring_roadmap/4cde_output/
"""

from __future__ import annotations

import gzip
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

from ntvetools import load_gtf_df


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
OUT_DIR = Path(__file__).parent / "4cde_output"
TEMP_DIR = OUT_DIR / "deseq2_temp"
QUANT_DIR = ROOT / "resources" / "quantfiles_filtered_pipeline" / "mus_musculus"
LIMS_FILE = first_existing_path(
    ROOT / "resources" / "Project_1716_lims_simplified.csv",
    ROOT / "Figure4" / "Project_1716_lims_simplified.csv",
)
GTF_FILE = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
DOCKER_DIR = ROOT / "resources" / "docker_DeSeq_with_Shrinkage"
DOCKER_IMAGE = "deseq2-with-shrinkage:latest"


def parse_sample_name(name: str) -> dict | None:
    match = re.match(r"S(\d+)_R(\d+)_TP(\d+)", str(name))
    if not match:
        return None
    return {
        "Sample": int(match.group(1)),
        "Replicate": int(match.group(2)),
        "TimePoint": int(match.group(3)),
    }


def load_gene_mappings() -> tuple[defaultdict, defaultdict, defaultdict]:
    print("Loading GTF reference data...")
    gtf_data = load_gtf_df(str(GTF_FILE))
    gtf_df = gtf_data["gtf_df"]
    gene_level_gtf = gtf_df.query("feature=='gene'")[["gene_id", "gene_biotype"]]
    gene_biotype = gene_level_gtf.set_index("gene_id")["gene_biotype"].to_dict()
    transcript_to_gene = gtf_df.set_index("transcript_id")["gene_id"].to_dict()
    gene_to_name = gtf_data["gene_id_to_name"]
    return (
        defaultdict(str, transcript_to_gene),
        defaultdict(str, gene_to_name),
        defaultdict(lambda: "NA", gene_biotype),
    )


def load_read_counts_for_sample(sample_file_path: Path, transcript_to_gene: dict) -> dict:
    with gzip.open(sample_file_path, "rt") as handle:
        df = pd.read_csv(handle, sep="\t")
    df["GeneID"] = df["Name"].str.split(".").str[0].map(lambda tid: transcript_to_gene[tid])
    return df.groupby("GeneID")["NumReads"].sum().to_dict()


def create_deseq2_r_script(sample_num: int) -> str:
    return f"""#!/usr/bin/env Rscript
library(DESeq2)

cat("\\n=== DESeq2: Sample {sample_num} - TP2 vs TP3 ===\\n")
cat("Method: Wald test with design ~ Condition and apeglm shrinkage\\n\\n")

setwd("/workspace/input")
output_dir <- "/workspace/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count_matrix <- read.csv("count_matrix.csv", row.names = 1, check.names = FALSE)
sample_metadata <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)
sample_metadata <- sample_metadata[order(match(sample_metadata$Sample, colnames(count_matrix))), ]
count_matrix <- count_matrix[, sample_metadata$Sample]

rownames(sample_metadata) <- sample_metadata$Sample
sample_metadata$Condition <- factor(sample_metadata$Condition, levels = c("Reference", "Target"))

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData = sample_metadata,
  design = ~ Condition
)

dds <- DESeq(dds, test="Wald")
res <- results(dds, contrast=c("Condition", "Target", "Reference"), alpha=0.05)

cat("Applying apeglm shrinkage...")
tryCatch({{
  res_shrunk <- lfcShrink(dds, coef="Condition_Target_vs_Reference", type="apeglm")
  cat(" OK\\n")
}}, error=function(e) {{
  cat(" FAILED (using unshrunken)\\n")
  res_shrunk <<- res
}})

gene_id_to_name_csv <- read.csv("gene_id_to_name.csv", row.names = 1)
gene_id_to_name <- setNames(gene_id_to_name_csv$GeneName, rownames(gene_id_to_name_csv))

results_df <- as.data.frame(res_shrunk)
results_df$Gene <- rownames(results_df)
results_df$GeneName <- gene_id_to_name[results_df$Gene]
results_df <- results_df[, c("Gene", "GeneName", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

write.csv(results_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)
write.csv(
  results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ],
  file.path(output_dir, "deseq2_results_significant.csv"),
  row.names = FALSE
)
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- data.frame(Gene = rownames(norm_counts), norm_counts)
write.csv(norm_counts_df, file.path(output_dir, "normalized_counts.csv"), row.names = FALSE)
"""


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_docker_image() -> None:
    inspect = subprocess.run(
        ["docker", "image", "inspect", DOCKER_IMAGE],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if inspect.returncode == 0:
        print(f"Using existing Docker image: {DOCKER_IMAGE}")
        return
    print(f"Building Docker image: {DOCKER_IMAGE}")
    run(["docker", "build", "-t", DOCKER_IMAGE, str(DOCKER_DIR)])


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    TEMP_DIR.mkdir(exist_ok=True)
    ensure_docker_image()

    transcript_to_gene, gene_to_name, gene_biotype = load_gene_mappings()

    print("Loading LIMS metadata...")
    lims_df = pd.read_csv(LIMS_FILE)
    sample_info = {}
    for _, row in lims_df.iterrows():
        parsed = parse_sample_name(row["Sample_NameLIMS"])
        if parsed:
            sample_info[row["Sample_Name"]] = parsed

    tp2_samples = [sid for sid, info in sample_info.items() if info["TimePoint"] == 2]
    tp3_samples = [sid for sid, info in sample_info.items() if info["TimePoint"] == 3]

    comparisons = {}
    for sample_num in [1, 2, 3, 4]:
        sample_tp2 = [s for s in tp2_samples if sample_info[s]["Sample"] == sample_num]
        sample_tp3 = [s for s in tp3_samples if sample_info[s]["Sample"] == sample_num]
        if sample_tp2 and sample_tp3:
            comparisons[sample_num] = {"tp2": sample_tp2, "tp3": sample_tp3, "all_samples": sample_tp2 + sample_tp3}

    print(f"Preparing count matrices for {len(comparisons)} sample groups...")

    sample_read_counts: dict[str, dict] = {}
    all_genes: set[str] = set()
    for comp in comparisons.values():
        for sample_name in comp["all_samples"]:
            sf_path = QUANT_DIR / f"{sample_name}.sf.gz"
            counts_dict = load_read_counts_for_sample(sf_path, transcript_to_gene)
            sample_read_counts[sample_name] = counts_dict
            all_genes.update(g for g in counts_dict if g)
    all_genes = sorted(all_genes)

    for sample_num, comp in comparisons.items():
        print(f"\nSample {sample_num}: building count matrix")
        sample_cols = comp["all_samples"]
        count_df = pd.DataFrame(
            {
                sample: [sample_read_counts[sample].get(gene_id, 0) for gene_id in all_genes]
                for sample in sample_cols
            },
            index=all_genes,
        ).astype(int)
        count_df = count_df[[gene_biotype[gid] == "protein_coding" for gid in count_df.index]]
        count_df = count_df[(count_df[sample_cols].sum(axis=1) > 0).values]

        metadata = []
        for sample_id in comp["tp2"]:
            metadata.append(
                {
                    "Sample": sample_id,
                    "Condition": "Reference",
                    "TimePoint": 2,
                    "Sample_Num": sample_info[sample_id]["Sample"],
                    "Replicate": sample_info[sample_id]["Replicate"],
                }
            )
        for sample_id in comp["tp3"]:
            metadata.append(
                {
                    "Sample": sample_id,
                    "Condition": "Target",
                    "TimePoint": 3,
                    "Sample_Num": sample_info[sample_id]["Sample"],
                    "Replicate": sample_info[sample_id]["Replicate"],
                }
            )
        metadata_df = pd.DataFrame(metadata)

        sample_temp_dir = TEMP_DIR / f"sample{sample_num}"
        sample_temp_dir.mkdir(exist_ok=True, parents=True)
        sample_output_dir = OUT_DIR / f"sample{sample_num}"
        sample_output_dir.mkdir(exist_ok=True)

        count_df.to_csv(sample_temp_dir / "count_matrix.csv")
        metadata_df.to_csv(sample_temp_dir / "sample_metadata.csv", index=False)
        pd.DataFrame({"GeneName": [gene_to_name.get(g, g) for g in count_df.index]}, index=count_df.index).to_csv(
            sample_temp_dir / "gene_id_to_name.csv"
        )
        (sample_temp_dir / "run_deseq2.R").write_text(create_deseq2_r_script(sample_num))

        run(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{sample_temp_dir.absolute()}:/workspace/input:rw",
                "-v",
                f"{sample_output_dir.absolute()}:/workspace/output:rw",
                DOCKER_IMAGE,
                "Rscript",
                "/workspace/input/run_deseq2.R",
            ]
        )

        results_path = sample_output_dir / "deseq2_results.csv"
        if not results_path.exists():
            raise FileNotFoundError(f"Expected results not found: {results_path}")
        print(f"  Saved: {results_path}")

    summary_rows = []
    for sample_num in sorted(comparisons):
        results_df = pd.read_csv(OUT_DIR / f"sample{sample_num}" / "deseq2_results.csv")
        n_total = int((results_df["padj"] < 0.05).sum())
        n_up = int(((results_df["padj"] < 0.05) & (results_df["log2FoldChange"] > 0)).sum())
        n_down = int(((results_df["padj"] < 0.05) & (results_df["log2FoldChange"] < 0)).sum())
        summary_rows.append(
            {"Sample": sample_num, "Total_Genes": len(results_df), "DEGs_padj_0.05": n_total, "Upregulated": n_up, "Downregulated": n_down}
        )
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "summary_statistics.csv", index=False)
    print(f"\nSaved: {OUT_DIR / 'summary_statistics.csv'}")
    print("Done.")


if __name__ == "__main__":
    main()
