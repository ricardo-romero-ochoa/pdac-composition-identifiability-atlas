#!/usr/bin/env Rscript
# Convert a PDAC scRNA/snRNA h5ad reference into the table format consumed by
# scripts/run_scrna_mapping.R. Usage:
#   Rscript scripts/convert_h5ad_to_tables.R input.h5ad data/external/scrna
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/convert_h5ad_to_tables.R input.h5ad [out_dir]", call. = FALSE)
infile <- args[[1]]
out_dir <- args[[2]] %||% "data/external/scrna"
if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("Install Bioconductor package zellkonverter.", call. = FALSE)
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Install Bioconductor package SummarizedExperiment.", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sce <- zellkonverter::readH5AD(infile)
assay <- if ("logcounts" %in% SummarizedExperiment::assayNames(sce)) "logcounts" else SummarizedExperiment::assayNames(sce)[1]
expr <- as.matrix(SummarizedExperiment::assay(sce, assay))
meta <- as.data.frame(SummarizedExperiment::colData(sce))
meta$cell <- colnames(expr)
readr::write_tsv(tibble::as_tibble(expr, rownames = "gene"), file.path(out_dir, "scrna_expression.tsv.gz"))
readr::write_tsv(meta, file.path(out_dir, "scrna_metadata.tsv"))
cat("Wrote ", file.path(out_dir, "scrna_expression.tsv.gz"), " and metadata.\n", sep = "")
