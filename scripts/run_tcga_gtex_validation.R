#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")
cfg <- read_config("config/atlas_config.yml")
print(run_tcga_gtex_validation(cfg))
