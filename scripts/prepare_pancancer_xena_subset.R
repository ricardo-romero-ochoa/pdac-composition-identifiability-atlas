#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(yaml)
})
source("R/00_utils.R")

cfg <- read_config("config/atlas_config.yml")
vcfg <- cfg$validation$tcga_gtex %||% list()
pcfg <- cfg$validation$pancancer_specificity %||% list()

expr_raw <- vcfg$xena_expression_file %||% "data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz"
pheno_raw <- vcfg$xena_phenotype_file %||% "data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz"
gene_list <- vcfg$requested_gene_file %||% "data/external/cache/tcga_gtex_requested_genes.txt"
probemap <- vcfg$xena_probemap_file %||% "data/external/tcga_gtex/gencode.v23.annotation.gene.probemap"
out_expr <- pcfg$expression_file %||% "data/external/pancancer/pancancer_expression.tsv.gz"
out_meta <- pcfg$metadata_file %||% "data/external/pancancer/pancancer_metadata.tsv"

missing <- c(expr_raw, pheno_raw, gene_list)[!file.exists(c(expr_raw, pheno_raw, gene_list))]
if (length(missing)) {
  stop(
    "Pan-cancer preparation reuses the local UCSC Xena Toil files from TCGA/GTEx validation. Missing: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

py <- Sys.which("python3")
if (!nzchar(py)) py <- Sys.which("python")
if (!nzchar(py)) stop("Python was not found on PATH.", call. = FALSE)

args <- c(
  "scripts/create_pancancer_xena_subset.py",
  "--expression", expr_raw,
  "--phenotype", pheno_raw,
  "--gene-list", gene_list,
  "--probemap", probemap,
  "--out-expression", out_expr,
  "--out-metadata", out_meta,
  "--out-counts", "data/external/pancancer/pancancer_cohort_counts.tsv"
)
cat("Building tissue-matched pan-cancer specificity subset from local Xena files...\n")
status <- system2(py, args)
if (!identical(status, 0L)) stop("Pan-cancer Xena subset builder failed with status ", status, call. = FALSE)
cat("Pan-cancer subset prepared. Now run:\n  Rscript scripts/run_pan_cancer_specificity_panel.R\n")
