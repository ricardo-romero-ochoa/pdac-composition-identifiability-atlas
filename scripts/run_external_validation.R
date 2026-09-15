#!/usr/bin/env Rscript
# Run the external-validation layer after the base atlas pipeline has finished.
# This script first exports the base atlas dataset_list target to an ordinary
# RDS cache file, then runs the optional external pipeline in its own targets
# store. Do not run the external script in the default _targets store.
# Equivalent low-level command after cache materialization:
#   Rscript -e "targets::tar_make(script = '_targets_external.R', store = '_targets_external')"

source("R/00_utils.R")
source("R/15_external_validation_utils.R")

cfg <- read_config("config/atlas_config.yml")
cat("Materializing base atlas target cache: dataset_list\n")
materialize_base_target_cache("dataset_list", cfg = cfg, store = "_targets", overwrite = TRUE)

cat("Running external validation pipeline in store: _targets_external\n")
targets::tar_make(script = "_targets_external.R", store = "_targets_external")
