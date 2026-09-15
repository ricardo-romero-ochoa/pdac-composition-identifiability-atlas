#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/01_download_geo.R")
source("R/02_metadata_curation.R")
cfg <- read_config()
esets <- download_all_geo(cfg)
export_metadata_audit(esets, cfg)
message("Curated metadata written to results/tables/curated_metadata_all.tsv")
