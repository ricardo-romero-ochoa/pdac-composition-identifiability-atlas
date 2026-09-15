#!/usr/bin/env Rscript
# Diagnose whether the external validation layer can see the GSE91035 transition signature.
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")

cfg <- read_config("config/atlas_config.yml")
catalog <- build_signature_catalog(cfg)
cat("Transition up genes:", length(catalog$transition_up$up), "
")
cat("Transition down genes:", length(catalog$transition_composite$down), "
")
cat("Transition source:", attr(catalog$tabs$transition, "source_file") %||% NA_character_, "
")

vcfg <- get_tcga_gtex_cfg(cfg)
expr_file <- vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"
if (file.exists(expr_file)) {
  expr <- read_table_auto(expr_file) |> standardize_gene_matrix()
  print(signature_catalog_diagnostics(expr, catalog))
} else {
  cat("TCGA/GTEx expression subset not found:", expr_file, "
")
}
