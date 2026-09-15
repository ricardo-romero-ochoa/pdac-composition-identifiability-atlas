#!/usr/bin/env Rscript
# Prepare compact scRNA/snRNA TSV files from an h5ad reference. This avoids
# reading h5ad from inside targets and avoids reticulate/zellkonverter Python
# environment ambiguity during run_scrna_mapping.R.

source("R/00_utils.R")
source("R/15_external_validation_utils.R")

cfg <- read_config("config/atlas_config.yml")
scfg <- get_validation_subcfg(cfg, "scrna")
h5ad <- scfg$h5ad_file %||% "data/external/scrna/scrna_reference.h5ad"
out_dir <- dirname(scfg$expression_file %||% "data/external/scrna/scrna_expression.tsv.gz")
gene_file <- file.path(out_dir, "scrna_gene_universe.tsv")
max_cells <- scfg$max_cells_per_celltype %||% 2000
seed <- cfg$project$random_seed %||% 20260531

if (!file.exists(h5ad)) {
  stop("h5ad file not found: ", h5ad, call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Building scRNA gene universe from atlas signature tables...\n")
catalog <- build_signature_catalog(cfg)
module_genes <- unique(unlist(lapply(catalog[setdiff(names(catalog), "tabs")], function(sig) {
  c(sig$up %||% character(), sig$down %||% character(), names(sig$weights %||% numeric()))
})))

hub_tbl <- load_atlas_signature_tables(cfg)$strict_hubs
hub_genes <- character()
gcol <- intersect(c("feature", "gene", "symbol"), colnames(hub_tbl))[1]
if (!is.na(gcol) && nrow(hub_tbl)) {
  hub_genes <- utils::head(normalise_gene_symbols(hub_tbl[[gcol]]), scfg$top_hub_genes %||% 40)
}

genes <- sort(unique(normalise_gene_symbols(c(module_genes, hub_genes))))
genes <- genes[nzchar(genes) & !is.na(genes)]
readr::write_tsv(tibble::tibble(gene = genes), gene_file)
cat("Wrote gene universe: ", gene_file, " (", length(genes), " genes)\n", sep = "")

python <- scfg$python %||% Sys.getenv("PDAC_PYTHON", unset = "")
if (!nzchar(python)) python <- Sys.which("python3")
if (!nzchar(python)) python <- Sys.which("python")
if (!nzchar(python)) stop("No Python executable found. Set PDAC_PYTHON or validation.scrna.python.", call. = FALSE)

cmd <- c(
  "scripts/convert_h5ad_to_tables.py",
  "--h5ad", h5ad,
  "--out-dir", out_dir,
  "--genes", gene_file,
  "--max-cells-per-celltype", as.character(max_cells),
  "--seed", as.character(seed)
)
if (!is.null(scfg$cell_type_column) && !is.na(scfg$cell_type_column) && nzchar(as.character(scfg$cell_type_column))) {
  cmd <- c(cmd, "--cell-type-column", as.character(scfg$cell_type_column))
}
if (!is.null(scfg$condition_column) && !is.na(scfg$condition_column) && nzchar(as.character(scfg$condition_column))) {
  cmd <- c(cmd, "--condition-column", as.character(scfg$condition_column))
}
if (!is.null(scfg$sample_column) && !is.na(scfg$sample_column) && nzchar(as.character(scfg$sample_column))) {
  cmd <- c(cmd, "--sample-column", as.character(scfg$sample_column))
}
if (isFALSE(scfg$fail_if_unknown_celltype %||% TRUE)) {
  cmd <- c(cmd, "--allow-unknown-celltype")
}
cat("Using Python: ", python, "\n", sep = "")
cat("Running: ", paste(shQuote(c(python, cmd)), collapse = " "), "\n", sep = "")
status <- system2(python, args = cmd)
if (!identical(status, 0L)) {
  stop(
    "h5ad conversion failed. Install required packages into this Python:\n  ",
    python, " -m pip install anndata h5py pandas scipy numpy\n",
    call. = FALSE
  )
}
cat("Done. Now run:\n  Rscript scripts/run_scrna_mapping.R\n", sep = "")
