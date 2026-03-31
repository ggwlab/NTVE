#!/usr/bin/env Rscript
# DESeq2 Analysis with Natural Cubic Splines: Cardio Development

library(DESeq2)
library(splines)
library(ggplot2)

cat("=== DESeq2 Cubic Splines Analysis: Cardio Development ===")
cat("\n\nQuestion: How does gene expression change over development time?\n")
cat("Method: Natural cubic spline basis functions (df=4)\n")
cat("Timepoints: Day 0-10 (11 timepoints, 3 replicates each)\n\n")

setwd("/workspace/input")
output_dir <- "/workspace/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Load data
cat("Step 1: Loading data...\n")
count_matrix <- read.csv("count_matrix.csv", row.names = 1, check.names = FALSE)
sample_metadata <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)

cat(paste("  Genes:", nrow(count_matrix), "\n"))
cat(paste("  Samples:", ncol(count_matrix), "\n"))
cat(paste("  Timepoints:", paste(unique(sample_metadata$Time), collapse=", "), "\n\n"))

# Step 2: Create spline basis
cat("Step 2: Creating natural cubic spline basis...\n")
# Natural cubic splines with specified degrees of freedom
spline_basis <- ns(sample_metadata$Time, df = 4)
colnames(spline_basis) <- paste0("spline", 1:ncol(spline_basis))

cat(paste("  Created", ncol(spline_basis), "spline basis functions\n\n"))

# Add splines to coldata
coldata <- cbind(sample_metadata, spline_basis)
rownames(coldata) <- coldata$Sample

# Convert Replicate to factor
coldata$Replicate <- factor(coldata$Replicate)

# Create design formula from spline columns + Replicate as fixed effect
spline_cols <- paste0("spline", 1:ncol(spline_basis))
design_formula <- as.formula(paste("~ Replicate +", paste(spline_cols, collapse=" + ")))
cat(paste("Design formula:", deparse(design_formula), "\n"))
cat("Note: Replicate included as fixed effect to account for replicate-specific variation\n\n")

# Step 3: Create DESeq2 object
cat("Step 3: Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData = coldata,
  design = design_formula
)

cat(paste("  Object created:", nrow(dds), "genes x", ncol(dds), "samples\n\n"))

# Step 4: Estimate size factors and dispersions
cat("Step 4: Estimating size factors and dispersions...\n")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds)

# Use less stringent outlier calling for time-course data
dds <- estimateDispersionsMAP(dds, outlierSD = 10)

cat(paste("  Dispersions estimated (outlierSD = 10 for time-course)\n\n"))

# Step 5: Perform LRT test
cat("Step 5: Performing Likelihood Ratio Test...\n")
cat(paste("  Full model:", deparse(design_formula), "\n"))
reduced_formula <- as.formula("~ Replicate")
cat(paste("  Reduced model:", deparse(reduced_formula), "\n"))
cat("  Testing: Is there significant time-dependent expression change (beyond replicate effects)?\n\n")

dds <- nbinomLRT(
  dds,
  full = design_formula,
  reduced = reduced_formula
)

# Step 6: Extract and examine results
cat("Step 6: Extracting results...\n")
res <- results(dds)

cat(paste("  Genes analyzed:", nrow(res), "\n"))
cat(paste("  Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm=TRUE), "\n"))
cat(paste("  Significant genes (padj < 0.1):", sum(res$padj < 0.1, na.rm=TRUE), "\n\n"))

# Show top genes
cat("Top 15 genes by adjusted p-value:\n")
top_genes <- res[order(res$padj), ][1:15, ]
print(top_genes[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")])
cat("\n")

# Step 7: Save results
cat("Step 7: Saving results...\n")

# Results table
results_df <- as.data.frame(res)
results_df$Gene <- rownames(results_df)
results_df <- results_df[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

results_file <- file.path(output_dir, "deseq2_results.csv")
write.csv(results_df, results_file, row.names = FALSE)
cat(paste("  Results saved:", results_file, "\n"))

# Significant genes only
sig_genes <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
sig_file <- file.path(output_dir, "deseq2_significant_genes.csv")
write.csv(sig_genes, sig_file, row.names = FALSE)
cat(paste("  Significant genes:", sig_file, "\n"))

# Save DESeq2 object for downstream analysis
dds_file <- file.path(output_dir, "dds.rds")
saveRDS(dds, dds_file)
cat(paste("  DESeq2 object:", dds_file, "\n"))

# Step 8: Create visualizations
cat("\nStep 8: Creating visualizations...\n")

# MA plot
pdf(file.path(output_dir, "ma_plot.pdf"), width = 8, height = 6)
plotMA(res, main="MA Plot: Log2 Fold Change vs Mean Count")
dev.off()
cat("  MA plot saved\n")

# P-value histogram
pdf(file.path(output_dir, "pvalue_histogram.pdf"), width = 8, height = 6)
hist(res$pvalue, breaks = 50, main = "P-value Distribution",
     xlab = "P-value", col = "lightblue", border = "white")
abline(v = 0.05, col = "red", lwd = 2, lty = 2)
dev.off()
cat("  P-value histogram saved\n")

# Dispersion plot
pdf(file.path(output_dir, "dispersion_plot.pdf"), width = 8, height = 6)
plotDispEsts(dds, main="Gene-wise Dispersion Estimates")
dev.off()
cat("  Dispersion plot saved\n")

cat("\n=== Analysis Complete ===")
cat("\nSummary:\n")
cat(paste("  Total genes analyzed:", nrow(res), "\n"))
cat(paste("  Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm=TRUE), "\n"))
cat(paste("  Mean baseMean expression:", round(mean(res$baseMean, na.rm=TRUE), 2), "\n"))
cat(paste("\nAll results saved to:", output_dir, "\n"))
