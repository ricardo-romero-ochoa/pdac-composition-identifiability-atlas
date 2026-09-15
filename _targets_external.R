library(targets)
library(tarchetypes)

source("R/00_utils.R")
source("R/04_de_limma.R")
source("R/07_deconvolution_adjustment.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")
source("R/17_scrna_mapping.R")
source("R/18_deconv_sensitivity.R")
source("R/19_depmap_opentargets.R")
source("R/20_external_validation_report.R")

tar_option_set(
  packages = c(
    "tidyverse", "data.table", "yaml", "matrixStats", "ggplot2", "readr",
    "targets", "survival", "AnnotationDbi", "org.Hs.eg.db", "httr2", "jsonlite"
  ),
  format = "rds"
)

list(
  tar_target(cfg, read_config("config/atlas_config.yml")),
  tar_target(external_tcga_gtex_validation, run_tcga_gtex_validation(cfg)),
  tar_target(external_scrna_mapping, run_scrna_signature_mapping(cfg)),
  tar_target(base_dataset_list_file, require_base_target_cache_file("dataset_list", cfg), format = "file"),
  tar_target(base_dataset_list, readRDS(base_dataset_list_file)),
  tar_target(external_deconv_sensitivity, run_deconv_sensitivity(base_dataset_list, cfg)),
  tar_target(external_target_prioritization, make_target_prioritization_table(cfg)),
  tar_target(external_validation_summary, run_external_validation_summary(cfg)),
  tar_target(external_validation_report, {
    external_tcga_gtex_validation
    external_scrna_mapping
    external_deconv_sensitivity
    external_target_prioritization
    external_validation_summary
    render_external_validation_report(quiet = TRUE)
  })
)
