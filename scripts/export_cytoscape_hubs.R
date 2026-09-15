#!/usr/bin/env Rscript
# Export the top hub-gene table as a Cytoscape-ready node table.

library(tidyverse)
source("R/00_utils.R")
cfg <- read_config()
hub_path <- file.path(cfg$project$output_dir, "tables", "integrative_hub_gene_prioritization.tsv")
if (!file.exists(hub_path)) stop("Run the main pipeline first.")

hubs <- readr::read_tsv(hub_path, show_col_types = FALSE) |>
  dplyr::slice_head(n = 200) |>
  dplyr::transmute(
    id = gene,
    label = gene,
    module = module,
    meta_logFC = meta_logFC,
    fdr = fdr,
    I2 = I2,
    integrated_hub_score = integrated_hub_score,
    core_flag = core_flag,
    microenvironment_flag = microenvironment_flag,
    transition_flag = transition_flag,
    retained_after_tme_adjustment = retained_after_tme_adjustment,
    tractability = tractability
  )
write_tsv(hubs, file.path(cfg$project$output_dir, "tables", "cytoscape_hub_nodes.tsv"))
message("Cytoscape node table written to results/tables/cytoscape_hub_nodes.tsv")
