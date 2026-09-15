#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")

cfg <- read_config("config/atlas_config.yml")
vcfg <- get_tcga_gtex_cfg(cfg)
cache_file <- vcfg$requested_gene_file %||% "data/external/cache/tcga_gtex_requested_genes.txt"
expr_file <- vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"
meta_file <- vcfg$metadata_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv"

files_to_remove <- unique(c(cache_file, expr_file, meta_file))
for (f in files_to_remove) {
  if (file.exists(f)) {
    cat("Removing stale file: ", f, "\n", sep = "")
    unlink(f)
  }
}

cat("Rebuilding requested gene universe and memory-safe TCGA/GTEx subset...\n")
status <- system2("Rscript", c("scripts/prepare_tcga_gtex_subset.R"))
if (!identical(status, 0L)) stop("prepare_tcga_gtex_subset.R failed with status ", status, call. = FALSE)

cat("Running cache diagnostic...\n")
status <- system2("Rscript", c("scripts/check_tcga_gtex_gene_cache.R"))
if (!identical(status, 0L)) stop("check_tcga_gtex_gene_cache.R failed with status ", status, call. = FALSE)

cat("Running external validation graph...\n")
status <- system2("Rscript", c("scripts/run_external_validation.R"))
if (!identical(status, 0L)) stop("run_external_validation.R failed with status ", status, call. = FALSE)

cat("Re-rendering external validation report if possible...\n")
status <- system2("Rscript", c("scripts/render_external_validation_report.R"))
if (!identical(status, 0L)) {
  warning("render_external_validation_report.R failed with status ", status, "; validation tables may still be present.")
}

cat("Done. Inspect results/tables/external/tcga_gtex_signature_gene_overlap.tsv and results/tables/external/tcga_survival_results.tsv.\n")
