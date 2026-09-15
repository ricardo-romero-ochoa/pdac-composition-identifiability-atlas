#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})
source("R/00_utils.R")
source("R/04_de_limma.R")
source("R/15_external_validation_utils.R")

cfg <- read_config("config/atlas_config.yml")
out_dir <- "results/external_deconvolution_inputs/cibersortx"
mix_dir <- file.path(out_dir, "mixtures")
dir.create(mix_dir, recursive = TRUE, showWarnings = FALSE)

bulk_cache <- "data/external/cache/base_dataset_list.rds"
if (!file.exists(bulk_cache)) {
  stop(
    "Missing cached base dataset list: ", bulk_cache,
    ". The base atlas has already been run in this project; do not redownload GEO. ",
    "Restore/copy the cache from the working repository.",
    call. = FALSE
  )
}

# Build/refesh the full-transcriptome 7-class single-cell reference.
h5ad <- cfg$validation$scrna$h5ad_file %||% "data/external/scrna/scrna_reference.h5ad"
py <- Sys.which("python3")
if (!nzchar(py)) py <- Sys.which("python")
if (!nzchar(py)) stop("Python was not found on PATH.", call. = FALSE)

ref_file <- file.path(out_dir, "PDAC_7class_refsample.txt")
if (!file.exists(ref_file) || is.na(file.info(ref_file)$size) || file.info(ref_file)$size == 0) {
  cat("Exporting full-transcriptome 7-class PDAC scRNA reference...\n")
  status <- system2(
    py,
    c(
      "scripts/export_cibersortx_pdac_reference.py",
      "--h5ad", h5ad,
      "--outdir", out_dir,
      "--max-cells-per-class", "700",
      "--min-expression-fraction", "0.01",
      "--max-genes", "18000",
      "--seed", "20260913"
    )
  )
  if (!identical(status, 0L)) {
    stop("CIBERSORTx scRNA-reference export failed with status ", status, call. = FALSE)
  }
} else {
  cat("Reusing existing CIBERSORTx reference: ", ref_file, "\n", sep = "")
}

if (!file.exists(ref_file) || file.info(ref_file)$size == 0) {
  stop("CIBERSORTx reference file was not created.", call. = FALSE)
}

datasets <- readRDS(bulk_cache)
if (!is.list(datasets) || !length(datasets)) {
  stop("Base dataset cache is empty or malformed.", call. = FALSE)
}

detect_log2 <- function(mat) {
  v <- as.numeric(mat)
  v <- v[is.finite(v)]
  if (!length(v)) return(FALSE)
  q <- stats::quantile(v, probs = c(0.01, 0.5, 0.99), na.rm = TRUE, names = FALSE)
  # RMA/normalized microarray matrices are typically on a compact log2 scale.
  isTRUE(q[3] < 35 && q[2] < 25)
}

row_median_impute <- function(mat) {
  n_missing <- sum(!is.finite(mat))
  if (!n_missing) return(list(mat = mat, n_imputed = 0L))
  for (i in seq_len(nrow(mat))) {
    bad <- !is.finite(mat[i, ])
    if (!any(bad)) next
    med <- suppressWarnings(stats::median(mat[i, !bad], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    mat[i, bad] <- med
  }
  list(mat = mat, n_imputed = n_missing)
}

manifest_rows <- list()
sample_rows <- list()

for (dsid in names(datasets)) {
  ds <- datasets[[dsid]]
  prepared <- average_technical_replicates(
    ds$expr,
    ds$meta,
    enabled = isTRUE(ds$config$technical_replicates)
  )
  expr <- standardize_gene_matrix(prepared$expr, collapse_fun = "mean")
  meta <- prepared$meta

  # CIBERSORTx mixture files must not contain missing/non-finite values.
  imp <- row_median_impute(expr)
  expr <- imp$mat

  was_log2 <- detect_log2(expr)
  transform <- "none_linear"
  if (was_log2) {
    expr <- 2^expr
    transform <- "inverse_log2_2powx"
  }

  if (any(!is.finite(expr))) {
    stop(dsid, ": non-finite values remain after mixture preparation.", call. = FALSE)
  }
  if (min(expr, na.rm = TRUE) < -1e-8) {
    stop(
      dsid, ": prepared mixture contains negative values after transformation. ",
      "CIBERSORTx requires a non-log/non-negative mixture matrix.",
      call. = FALSE
    )
  }

  # Deterministic gene ordering.
  ord <- order(rownames(expr))
  expr <- expr[ord, , drop = FALSE]
  out <- as.data.frame(expr, check.names = FALSE)
  out <- cbind(Gene = rownames(expr), out)
  out_file <- file.path(mix_dir, paste0(dsid, "_mixture.txt"))
  write_tsv(out, out_file)

  manifest_rows[[dsid]] <- tibble(
    dataset = dsid,
    mixture_file = out_file,
    n_genes = nrow(expr),
    n_samples = ncol(expr),
    transformation = transform,
    n_values_imputed = imp$n_imputed,
    min_value = min(expr, na.rm = TRUE),
    median_value = stats::median(expr, na.rm = TRUE),
    max_value = max(expr, na.rm = TRUE)
  )

  if (!is.null(meta) && nrow(meta)) {
    sample_col <- intersect(c("sample", "sample_id", "array", "gsm"), names(meta))[1]
    condition_col <- intersect(c("condition", "group", "class", "phenotype"), names(meta))[1]
    if (!is.na(sample_col)) {
      x <- tibble(
        dataset = dsid,
        sample = as.character(meta[[sample_col]])
      )
      if (!is.na(condition_col)) x$condition <- as.character(meta[[condition_col]])
      sample_rows[[dsid]] <- x
    }
  }

  cat(
    dsid, ": wrote ", nrow(expr), " genes x ", ncol(expr),
    " samples; transform=", transform,
    "; imputed=", imp$n_imputed, "\n",
    sep = ""
  )
}

manifest <- bind_rows(manifest_rows)
write_tsv(manifest, file.path(out_dir, "PDAC_CIBERSORTx_mixture_manifest.tsv"))
if (length(sample_rows)) {
  write_tsv(bind_rows(sample_rows), file.path(out_dir, "PDAC_CIBERSORTx_sample_manifest.tsv"))
}

required <- c(
  "GSE15471", "GSE28735", "GSE62165",
  "GSE16515", "GSE71989", "GSE91035"
)
missing <- setdiff(required, manifest$dataset)
if (length(missing)) {
  stop("Missing required PDAC cohorts from CIBERSORTx export: ", paste(missing, collapse = ", "), call. = FALSE)
}

cat("\nCIBERSORTx final-validation input package is ready:\n  ", out_dir, "\n", sep = "")
cat("Reference: ", ref_file, "\n", sep = "")
cat("Mixtures:  ", length(required), " cohort files\n", sep = "")
cat("\nNext external step:\n  bash scripts/run_cibersortx_final_validation.sh\n")
