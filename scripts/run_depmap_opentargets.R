#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/19_depmap_opentargets.R")
cfg <- read_config("config/atlas_config.yml")
print(make_target_prioritization_table(cfg))
