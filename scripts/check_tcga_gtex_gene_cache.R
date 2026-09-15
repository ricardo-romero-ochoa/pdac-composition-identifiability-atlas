#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")

args <- commandArgs(trailingOnly = TRUE)
fail_on_mismatch <- !any(args %in% c("--no-fail", "--warn-only"))

cfg <- read_config("config/atlas_config.yml")
vcfg <- get_tcga_gtex_cfg(cfg)
out_dir <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_file <- vcfg$requested_gene_file %||% "data/external/cache/tcga_gtex_requested_genes.txt"
expr_file <- vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"

cat("Checking TCGA/GTEx requested-gene cache\n")
cat("Cache file: ", cache_file, "\n", sep = "")
cat("Subset matrix: ", expr_file, "\n", sep = "")

catalog <- build_signature_catalog(cfg)
expected <- sort(unique(signature_gene_universe(catalog)))
expected <- expected[!is.na(expected) & nzchar(expected)]

cached <- character()
if (file.exists(cache_file)) {
  cached <- readLines(cache_file, warn = FALSE)
  cached <- sort(unique(normalise_gene_symbols(cached)))
  cached <- cached[!is.na(cached) & nzchar(cached)]
}

missing_from_cache <- setdiff(expected, cached)
extra_in_cache <- setdiff(cached, expected)

# IMPORTANT: compute the match flag before constructing the tibble.
# tibble() uses data masking and evaluates columns sequentially; if this
# expression is placed inside tibble() after columns named missing_from_cache
# and extra_in_cache are created, those names resolve to the scalar count
# columns (length 1) instead of the vectors above. That makes a perfect cache
# match (0 missing, 0 extra) evaluate FALSE.
n_missing_from_cache <- length(missing_from_cache)
n_extra_in_cache <- length(extra_in_cache)
cache_matches_signature_universe <-
  n_missing_from_cache == 0L &&
  n_extra_in_cache == 0L &&
  length(expected) > 0L

subset_genes <- character()
if (file.exists(expr_file)) {
  con <- if (grepl("\\.gz$", expr_file)) gzfile(expr_file, "rt") else file(expr_file, "rt")
  on.exit(close(con), add = TRUE)
  hdr <- readLines(con, n = 1, warn = FALSE)
  # The subset matrix is expected to have genes as rows in the first column.
  # Full parsing can be expensive, so this diagnostic only reports whether the
  # file exists. Coverage is validated by run_tcga_gtex_validation.R.
}

summary <- tibble(
  expected_signature_genes = length(expected),
  cached_genes = length(cached),
  missing_from_cache = n_missing_from_cache,
  extra_in_cache = n_extra_in_cache,
  cache_file_exists = file.exists(cache_file),
  subset_matrix_exists = file.exists(expr_file),
  cache_matches_signature_universe = cache_matches_signature_universe
)
write_tsv(summary, file.path(out_dir, "tcga_gtex_gene_cache_audit.tsv"))
writeLines(missing_from_cache, file.path(out_dir, "tcga_gtex_genes_missing_from_cache.txt"), useBytes = TRUE)
writeLines(extra_in_cache, file.path(out_dir, "tcga_gtex_genes_extra_in_cache.txt"), useBytes = TRUE)

print(summary)
cat("Wrote diagnostics to ", out_dir, "\n", sep = "")

if (fail_on_mismatch && !isTRUE(cache_matches_signature_universe)) {
  stop("TCGA/GTEx requested-gene cache does not match the current signature universe. Run scripts/refresh_tcga_gtex_cache_and_validation.R", call. = FALSE)
}
