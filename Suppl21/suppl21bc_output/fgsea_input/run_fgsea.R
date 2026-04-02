#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(BiocParallel)
  library(jsonlite)
})

register(MulticoreParam(workers = max(1, parallel::detectCores() - 1)))

setwd("/workspace/input")
output_dir <- "/workspace/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_SIZE <- 15
MAX_SIZE <- 500

c2_json <- fromJSON("/workspace/input/c2.all.v2025.1.Hs.json")
c2_pathways <- lapply(names(c2_json), function(name) c2_json[[name]]$geneSymbols)
names(c2_pathways) <- names(c2_json)

c5_json <- fromJSON("/workspace/input/c5.go.v2025.1.Hs.json")
c5_pathways <- lapply(names(c5_json), function(name) c5_json[[name]]$geneSymbols)
names(c5_pathways) <- names(c5_json)

run_fgsea_timepoint <- function(rnk_file, pathways) {
  timepoint <- gsub("ranked_genes_day(\\d+)\\.rnk", "\\1", basename(rnk_file))
  ranks_df <- read.table(rnk_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  ranks <- setNames(ranks_df$rank, ranks_df$gene)
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = ranks, minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0.0)
  res <- as.data.frame(res)
  res <- res[order(res$pval), ]
  res$timepoint <- as.integer(timepoint)
  res
}

rnk_files <- sort(list.files(pattern="^ranked_genes_day.*\\.rnk$", full.names=TRUE))
c2_results_list <- lapply(rnk_files, function(f) run_fgsea_timepoint(f, c2_pathways))
c5_results_list <- lapply(rnk_files, function(f) run_fgsea_timepoint(f, c5_pathways))
c2_all <- do.call(rbind, c2_results_list)
c5_all <- do.call(rbind, c5_results_list)
fwrite(c2_all, file.path(output_dir, "fgsea_c2_all_timepoints.csv"), sep="\t")
fwrite(c5_all, file.path(output_dir, "fgsea_c5_all_timepoints.csv"), sep="\t")
