#!/usr/bin/env Rscript
# Exports gene-symbol expression matrices suitable for upload to CIBERSORT/CIBERSORTx.
# This script intentionally does not run CIBERSORT because LM22/CIBERSORTx usage may
# require external files, licensing, or web/API credentials.

source("R/00_utils.R")
source("R/01_download_geo.R")
source("R/02_metadata_curation.R")
source("R/03_preprocess_microarray.R")

cfg <- read_config()
esets <- download_all_geo(cfg)
metadata <- curate_all_metadata(esets, cfg)
datasets <- preprocess_all(esets, metadata, cfg)

outdir <- file.path(cfg$project$output_dir, "cibersort_inputs")
ensure_dir(outdir)
for (dsid in names(datasets)) {
  expr <- datasets[[dsid]]$expr
  out <- data.frame(GeneSymbol = rownames(expr), expr, check.names = FALSE)
  data.table::fwrite(out, file.path(outdir, paste0(dsid, "_cibersort_input.tsv")), sep = "\t")
}
message("CIBERSORT-compatible inputs written to ", outdir)
