#!/usr/bin/env Rscript
library(DESeq2)

cat("\n=== DESeq2: Sample 2 - TP2 vs TP3 ===\n")
cat("Method: Wald test with design ~ Condition and apeglm shrinkage\n\n")

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
tryCatch({
  res_shrunk <- lfcShrink(dds, coef="Condition_Target_vs_Reference", type="apeglm")
  cat(" OK\n")
}, error=function(e) {
  cat(" FAILED (using unshrunken)\n")
  res_shrunk <<- res
})

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
