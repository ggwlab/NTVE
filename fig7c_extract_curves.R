#!/usr/bin/env Rscript
# Extract fitted smooth curves from DESeq2 cubic spline model: Cardio Development
#
# Run inside the repo's Docker container via fig7c_run.sh — do not run directly.
# Container paths:
#   /workspace/output  → Figure7/cardio_deseq2_cubicsplines_output/
#
# Requires: /workspace/output/dds.rds  (produced by run_deseq2_analysis.R)
# Outputs:  /workspace/output/fitted_trajectories_all_genes.csv
#           /workspace/output/glm_coefficients.csv
#           /workspace/output/fitted_curves_comparison.pdf
#           /workspace/output/top_genes_summary.csv

library(DESeq2)

cat("=== Extracting Fitted Curves from DESeq2 Spline Model (Cardio) ===\n\n")

output_dir <- "/workspace/output"

dds <- readRDS(file.path(output_dir, "dds.rds"))
cat(sprintf("Loaded DESeq2 object: %d genes x %d samples\n\n", nrow(dds), ncol(dds)))

# Step 1: GLM coefficients
cat("Step 1: Extracting GLM coefficients...\n")
coefs    <- coef(dds)
coefs_df <- as.data.frame(coefs)
coefs_df$Gene <- rownames(coefs_df)
write.csv(coefs_df, file.path(output_dir, "glm_coefficients.csv"), row.names = FALSE)
cat(sprintf("  %d genes x %d coefficients\n\n", nrow(coefs), ncol(coefs)))

# Step 2: Design matrix
cat("Step 2: Extracting design matrix...\n")
design_matrix <- model.matrix(design(dds), colData(dds))
cat(sprintf("  %d samples x %d columns\n\n", nrow(design_matrix), ncol(design_matrix)))

# Step 3: Metadata and normalised counts
cat("Step 3: Loading metadata and normalised counts...\n")
metadata    <- as.data.frame(colData(dds))
timepoints  <- sort(unique(metadata$Time))
norm_counts <- counts(dds, normalized = TRUE)
cat(sprintf("  Timepoints: %s\n\n", paste(timepoints, collapse = ", ")))

# Step 4: Fitted values (DESeq2 coefficients are on log2 scale)
cat("Step 4: Computing fitted values...\n")
fitted_values <- 2^(design_matrix %*% t(coefs))
cat(sprintf("  %d samples x %d genes\n\n", nrow(fitted_values), ncol(fitted_values)))

# Step 5: Load statistical results
cat("Step 5: Loading statistical results...\n")
results_file <- file.path(output_dir, "deseq2_results.csv")
if (!file.exists(results_file)) stop("deseq2_results.csv not found — run run_deseq2_analysis.R first.")
results   <- read.csv(results_file)
results   <- results[order(results$padj), ]
top_genes <- head(results[!is.na(results$padj), ], 20)
cat(sprintf("  Total genes: %d\n\n", nrow(results)))

# Step 6: Extract trajectories for ALL genes
cat("Step 6: Computing mean trajectories for all genes...\n")
genes <- rownames(coefs)
all_trajectories <- vector("list", length(genes) * length(timepoints))
idx <- 1L

for (i in seq_along(genes)) {
  gene          <- genes[i]
  gene_fitted   <- fitted_values[, gene]
  gene_observed <- norm_counts[gene, ]
  gene_padj_val <- results$padj[match(gene, results$Gene)]

  for (tp in timepoints) {
    rows <- metadata$Time == tp
    all_trajectories[[idx]] <- data.frame(
      Gene        = gene,
      Time        = tp,
      Fitted      = mean(gene_fitted[rows]),
      Fitted_SD   = sd(gene_fitted[rows]),
      Observed    = mean(gene_observed[rows]),
      Observed_SD = sd(gene_observed[rows]),
      padj        = gene_padj_val
    )
    idx <- idx + 1L
  }

  if (i %% 2000 == 0) cat(sprintf("  Processed %d / %d genes\n", i, length(genes)))
}

all_traj_df <- do.call(rbind, all_trajectories)
cat(sprintf("  Done: %d genes, %d timepoints\n\n", length(genes), length(timepoints)))

# Step 7: Save trajectories
cat("Step 7: Saving trajectory data...\n")
traj_file <- file.path(output_dir, "fitted_trajectories_all_genes.csv")
write.csv(all_traj_df, traj_file, row.names = FALSE)
cat(sprintf("  %s (%d rows)\n\n", traj_file, nrow(all_traj_df)))

# Step 8: Diagnostic plot (top 20 genes)
cat("Step 8: Diagnostic plot for top 20 genes...\n")
top_traj <- all_traj_df[all_traj_df$Gene %in% top_genes$Gene, ]

pdf(file.path(output_dir, "fitted_curves_comparison.pdf"), width = 16, height = 20)
par(mfrow = c(5, 4), mar = c(4, 4, 3, 1))
for (gene in top_genes$Gene) {
  gd <- top_traj[top_traj$Gene == gene, ]
  gd <- gd[order(gd$Time), ]
  y_range <- range(c(gd$Observed + gd$Observed_SD,
                     gd$Observed - gd$Observed_SD,
                     gd$Fitted), na.rm = TRUE)
  plot(gd$Time, gd$Observed, type = "o", pch = 16, col = "blue", lwd = 2,
       ylim = y_range, xlab = "Time (days)", ylab = "Normalised counts",
       main = sprintf("%s (padj=%.2e)", gene, gd$padj[1]))
  arrows(gd$Time, gd$Observed - gd$Observed_SD,
         gd$Time, gd$Observed + gd$Observed_SD,
         angle = 90, code = 3, length = 0.05, col = "blue")
  lines(gd$Time, gd$Fitted, type = "o", pch = 15, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Observed", "Fitted"), col = c("blue", "red"),
         lty = c(1, 2), pch = c(16, 15), lwd = 2, cex = 0.8)
  grid()
}
dev.off()

# Step 9: Summary table
top_genes_coefs <- merge(top_genes, coefs_df, by = "Gene")
write.csv(top_genes_coefs, file.path(output_dir, "top_genes_summary.csv"), row.names = FALSE)

cat("=== Complete ===\n")
cat("Next: run fig7c.py to produce the heatmap figures.\n")
