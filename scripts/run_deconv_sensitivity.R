#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/04_de_limma.R")
source("R/07_deconvolution_adjustment.R")
source("R/15_external_validation_utils.R")
source("R/18_deconv_sensitivity.R")

cfg <- read_config("config/atlas_config.yml")
if (!requireNamespace("targets", quietly = TRUE)) stop("targets is required. Run install_packages.R first.")
materialize_base_target_cache("dataset_list", cfg = cfg, store = "_targets", overwrite = TRUE)
dataset_list <- read_base_target_cache("dataset_list", cfg = cfg)
print(run_deconv_sensitivity(dataset_list, cfg))
