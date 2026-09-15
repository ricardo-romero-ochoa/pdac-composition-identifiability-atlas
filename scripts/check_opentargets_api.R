#!/usr/bin/env Rscript
# Smoke test for the Open Targets GraphQL layer.
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/19_depmap_opentargets.R")

cfg <- read_config("config/atlas_config.yml")
gene <- commandArgs(trailingOnly = TRUE)
gene <- if (length(gene)) gene[[1]] else "SLC6A14"
ann <- symbol_to_ensembl(gene)
print(ann)
ens <- ann$ensembl[[1]]
if (is.na(ens) || !nzchar(ens)) stop("Could not map symbol to Ensembl. Check org.Hs.eg.db.", call. = FALSE)
dat <- opentargets_query_target(ens, timeout_sec = get_target_cfg(cfg)$opentargets_timeout_sec %||% 30)
res <- summarise_opentargets_response(gene, ens, dat, disease_id = get_target_cfg(cfg)$opentargets_disease_id %||% NULL)
print(res)
cat("\nPancreas-specific disease score:", res$pancreas_disease_score %||% NA, "\n")
cat("Pancreas-specific known drug count:", res$pancreas_known_drug_count %||% NA, "\n")
if (!nrow(res) || !identical(res$opentargets_status[[1]], "ok")) quit(status = 1)
