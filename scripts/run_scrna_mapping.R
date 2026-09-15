#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/17_scrna_mapping.R")
cfg <- read_config("config/atlas_config.yml")
print(run_scrna_signature_mapping(cfg))
