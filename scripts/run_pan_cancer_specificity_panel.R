#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})
source("R/00_utils.R")
source("R/15_external_validation_utils.R")

rank_auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]; label <- as.integer(label[ok])
  n1 <- sum(label == 1); n0 <- sum(label == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cfg <- read_config("config/atlas_config.yml")
pcfg <- cfg$validation$pancancer_specificity %||% list()
out_dir <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- pcfg$expression_file %||% "data/external/pancancer/pancancer_expression.tsv.gz"
meta_file <- pcfg$metadata_file %||% "data/external/pancancer/pancancer_metadata.tsv"
if (!file.exists(expr_file) || !file.exists(meta_file)) {
  dir.create(dirname(expr_file), recursive = TRUE, showWarnings = FALSE)
  write_tsv(tibble(gene = c("TP53", "KRAS"), sample_1 = c(1, 2), sample_2 = c(3, 4)), sub("\\.gz$", "", expr_file))
  write_tsv(tibble(sample = c("sample_1", "sample_2"), cancer_type = c("COAD", "COAD"), class = c("tumor", "gtex_normal")), meta_file)
  stop("Pan-cancer expression/metadata files are missing. Template files were written under data/external/pancancer/. Provide a gene-by-sample matrix and sample metadata, then rerun.", call. = FALSE)
}

cat("Reading pan-cancer expression and metadata...\n")
expr <- read_tsv(expr_file, show_col_types = FALSE)
gene_col <- intersect(c("gene", "symbol", "feature"), names(expr))[1]
if (is.na(gene_col)) stop("Expression matrix must have a gene/symbol/feature column.", call. = FALSE)
mat_df <- as.data.frame(expr, check.names = FALSE)
genes <- normalise_gene_symbols(mat_df[[gene_col]])
mat_df[[gene_col]] <- NULL
mat <- as.matrix(as.data.frame(lapply(mat_df, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE))
rownames(mat) <- genes
keep <- !is.na(rownames(mat)) & nzchar(rownames(mat))
mat <- mat[keep, , drop = FALSE]
meta <- read_tsv(meta_file, show_col_types = FALSE)
sample_col <- pcfg$sample_column %||% "sample"
type_col <- pcfg$cancer_type_column %||% "cancer_type"
class_col <- pcfg$class_column %||% "class"
missing_cols <- setdiff(c(sample_col, type_col, class_col), names(meta))
if (length(missing_cols)) stop("Pan-cancer metadata missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
common <- intersect(colnames(mat), meta[[sample_col]])
if (!length(common)) stop("No overlapping samples between pan-cancer matrix and metadata.", call. = FALSE)
mat <- mat[, common, drop = FALSE]
meta <- meta[match(common, meta[[sample_col]]), , drop = FALSE]

catalog <- build_signature_catalog(cfg)
scores <- score_catalog_on_matrix(mat, catalog, min_genes = 3)
mods <- setdiff(colnames(scores), "sample")
meta_scores <- left_join(meta, scores, by = setNames("sample", sample_col))

tumor_label <- pcfg$tumor_label %||% "tumor"
normal_label <- pcfg$normal_label %||% "gtex_normal"
min_t <- pcfg$min_tumor %||% 20
min_n <- pcfg$min_normal %||% 20

res <- meta_scores |>
  filter(.data[[class_col]] %in% c(tumor_label, normal_label)) |>
  mutate(cancer_type_for_panel = .data[[type_col]]) |>
  group_by(.data$cancer_type_for_panel) |>
  group_modify(function(df, key) {
    n_t <- sum(df[[class_col]] == tumor_label)
    n_n <- sum(df[[class_col]] == normal_label)
    if (n_t < min_t || n_n < min_n) return(tibble())
    bind_rows(lapply(mods, function(m) {
      s <- df[[m]]
      lab <- df[[class_col]] == tumor_label
      auc <- rank_auc(s, lab)
      tibble(module = m, n_tumor = n_t, n_normal = n_n, auc = auc,
             reversed_auc = ifelse(is.na(auc), NA_real_, max(auc, 1 - auc)),
             pdac_specificity_flag = ifelse(!is.na(auc) & auc >= 0.90, "fails_specificity_if_non_pdac", "ok"))
    }))
  }) |>
  ungroup() |>
  rename(cancer_type = cancer_type_for_panel)

write_tsv(res, file.path(out_dir, "pancancer_specificity_panel.tsv"))
cat("Wrote ", nrow(res), " pan-cancer specificity rows.\n", sep = "")
