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

run_fgsea_timepoint <- function(rnk_file, pathways, pathway_name, min_size, max_size) {
  ranks_df <- fread(rnk_file)
  ranks <- setNames(ranks_df$rank, ranks_df$gene)
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  res <- fgsea(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size)
  tp <- sub(".*day([0-9]+)\\.rnk$", "\\1", rnk_file)
  res$timepoint <- as.integer(tp)
  res$collection <- pathway_name
  res
}

rnk_files <- list.files("/workspace/input", pattern = "^ranked_genes_day[0-9]+\\.rnk$", full.names = TRUE)
rnk_files <- sort(rnk_files)

c2_res <- rbindlist(lapply(rnk_files, run_fgsea_timepoint, pathways = c2_pathways, pathway_name = "C2", min_size = MIN_SIZE, max_size = MAX_SIZE), fill = TRUE)
c5_res <- rbindlist(lapply(rnk_files, run_fgsea_timepoint, pathways = c5_pathways, pathway_name = "C5", min_size = MIN_SIZE, max_size = MAX_SIZE), fill = TRUE)

fwrite(c2_res, file.path(output_dir, "fgsea_c2_all_timepoints.csv"), sep = "\t")
fwrite(c5_res, file.path(output_dir, "fgsea_c5_all_timepoints.csv"), sep = "\t")
