#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")

cfg <- read_config("config/atlas_config.yml")
vcfg <- get_tcga_gtex_cfg(cfg)
out <- vcfg$requested_gene_file %||% "data/external/cache/tcga_gtex_requested_genes.txt"

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
cat("Building atlas signature gene universe for TCGA/GTEx subsetting...\n")
catalog <- build_signature_catalog(cfg)
genes <- signature_gene_universe(catalog)
genes <- sort(unique(genes[!is.na(genes) & nzchar(genes)]))
writeLines(genes, out, useBytes = TRUE)
cat("Wrote ", length(genes), " genes to ", out, "\n", sep = "")
