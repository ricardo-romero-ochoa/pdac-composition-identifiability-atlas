#!/usr/bin/env Rscript
# Sample-level PDAC subtype benchmark using the published PurIST single-sample
# classifier (Rashid et al., Clin Cancer Res 2020).
#
# The published model is applied without retraining:
# score = -6.815 + sum_i beta_i I(GeneA_i > GeneB_i)
# probability_basal = plogis(score); >0.5 = basal-like.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})
source("R/00_utils.R")
source("R/15_external_validation_utils.R")

rank_auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]
  label <- as.integer(label[ok])
  n1 <- sum(label == 1L)
  n0 <- sum(label == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cfg <- read_config("config/atlas_config.yml")
out_dir <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_file <- file.path(out_dir, "tcga_gtex_module_scores.tsv")
if (!file.exists(module_file)) {
  stop("Missing PurIST benchmark module-score input: ", module_file, call. = FALSE)
}

purist_pairs <- tribble(
  ~pair, ~gene_a,    ~gene_b,    ~coefficient,
   1L,   "GPR87",    "REG4",       1.994,
   2L,   "KRT6A",    "ANXA10",     2.031,
   3L,   "BCAR3",    "GATA6",      1.618,
   4L,   "PTGES",    "CLDN18",     0.922,
   5L,   "ITGA3",    "LGALS4",     1.059,
   6L,   "C16ORF74", "DDC",        0.929,
   7L,   "S100A2",   "SLC40A1",    2.505,
   8L,   "KRT5",     "CLRN3",      0.485
)
purist_intercept <- -6.815
needed <- sort(unique(c(purist_pairs$gene_a, purist_pairs$gene_b)))

atlas_expr_file <- cfg$validation$tcga_gtex$expression_file %||%
  "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"
atlas_meta_file <- cfg$validation$tcga_gtex$metadata_file %||%
  "data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv"

purist_gene_file <- "data/external/cache/purist_requested_genes.txt"
purist_expr_file <- "data/external/tcga_gtex/purist_tcga_paad_expression.tsv.gz"
purist_meta_file <- "data/external/tcga_gtex/purist_tcga_paad_metadata.tsv"

dir.create(dirname(purist_gene_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(purist_expr_file), recursive = TRUE, showWarnings = FALSE)
writeLines(needed, purist_gene_file)

subset_has_all_genes <- function(path, genes) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
    return(FALSE)
  }
  ok <- tryCatch({
    x <- standardize_gene_matrix(
      read_tsv(path, show_col_types = FALSE),
      collapse_fun = "mean"
    )
    all(genes %in% rownames(x))
  }, error = function(e) FALSE)
  isTRUE(ok)
}

# Prefer the existing atlas subset only if it happens to contain all PurIST
# genes. Normally it will not, because the atlas cache is intentionally limited
# to the atlas signature universe.
use_atlas_subset <- subset_has_all_genes(atlas_expr_file, needed) &&
  file.exists(atlas_meta_file)

if (use_atlas_subset) {
  expr_file <- atlas_expr_file
  meta_file <- atlas_meta_file
  input_source <- "atlas_tcga_gtex_subset"
} else {
  if (!subset_has_all_genes(purist_expr_file, needed) || !file.exists(purist_meta_file)) {
    raw_expr <- cfg$validation$tcga_gtex$xena_expression_file %||%
      "data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz"
    raw_pheno <- cfg$validation$tcga_gtex$xena_phenotype_file %||%
      "data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz"
    probemap <- cfg$validation$tcga_gtex$xena_probemap_file %||%
      "data/external/tcga_gtex/gencode.v23.annotation.gene.probemap"

    missing_raw <- c(raw_expr, raw_pheno, probemap)[
      !file.exists(c(raw_expr, raw_pheno, probemap))
    ]
    if (length(missing_raw)) {
      stop(
        "The atlas TCGA/GTEx subset does not contain all 16 PurIST genes, and ",
        "the local raw Xena files needed to build the dedicated PurIST subset ",
        "are missing: ", paste(missing_raw, collapse = ", "),
        call. = FALSE
      )
    }

    py <- Sys.which("python3")
    if (!nzchar(py)) py <- Sys.which("python")
    if (!nzchar(py)) stop("Python was not found on PATH.", call. = FALSE)

    cat(
      "Atlas TCGA/GTEx subset is signature-restricted; building dedicated ",
      "16-gene PurIST subset from local Xena data...\n",
      sep = ""
    )

    args <- c(
      "scripts/create_tcga_gtex_xena_subset.py",
      "--expression", raw_expr,
      "--phenotype", raw_pheno,
      "--gene-list", purist_gene_file,
      "--probemap", probemap,
      "--out-expression", purist_expr_file,
      "--out-metadata", purist_meta_file,
      "--max-samples",
      as.character(cfg$validation$tcga_gtex$max_xena_samples %||% 1000L),
      "--no-download-probemap"
    )
    status <- system2(py, args)
    if (!identical(status, 0L)) {
      stop("Dedicated PurIST Xena subset builder failed with status ", status, call. = FALSE)
    }
  }

  expr_file <- purist_expr_file
  meta_file <- purist_meta_file
  input_source <- "dedicated_16_gene_xena_subset"
}

expr <- standardize_gene_matrix(
  read_tsv(expr_file, show_col_types = FALSE),
  collapse_fun = "mean"
)
meta <- read_tsv(meta_file, show_col_types = FALSE)
scores <- read_tsv(module_file, show_col_types = FALSE)

missing_genes <- setdiff(needed, rownames(expr))
pair_audit <- purist_pairs |>
  mutate(
    gene_a_available = .data$gene_a %in% rownames(expr),
    gene_b_available = .data$gene_b %in% rownames(expr),
    pair_available = .data$gene_a_available & .data$gene_b_available,
    expression_source = input_source
  )
write_tsv(pair_audit, file.path(out_dir, "purist_gene_pair_audit.tsv"))

if (length(missing_genes)) {
  stop(
    "PurIST dedicated input is still missing genes after Xena extraction: ",
    paste(missing_genes, collapse = ", "),
    call. = FALSE
  )
}

if (!all(c("sample", "condition") %in% names(meta))) {
  stop("TCGA/GTEx metadata must contain sample and condition.", call. = FALSE)
}
tumor_samples <- meta |>
  filter(tolower(as.character(.data$condition)) == "tumor") |>
  pull(.data$sample) |>
  as.character()
tumor_samples <- intersect(tumor_samples, colnames(expr))
if (length(tumor_samples) < 50L) {
  stop("Too few TCGA-PAAD tumors for PurIST: ", length(tumor_samples), call. = FALSE)
}

pair_mat <- sapply(seq_len(nrow(purist_pairs)), function(i) {
  a <- purist_pairs$gene_a[[i]]
  b <- purist_pairs$gene_b[[i]]
  as.numeric(expr[a, tumor_samples, drop = TRUE] > expr[b, tumor_samples, drop = TRUE])
})
pair_mat <- as.matrix(pair_mat)
rownames(pair_mat) <- tumor_samples
colnames(pair_mat) <- paste0("TSP", purist_pairs$pair)

linear_score <- purist_intercept + as.numeric(pair_mat %*% purist_pairs$coefficient)
prob_basal <- stats::plogis(linear_score)
subtype <- ifelse(prob_basal > 0.5, "basal_like", "classical")
graded <- cut(
  prob_basal,
  breaks = c(-Inf, 0.1, 0.5, 0.9, Inf),
  labels = c("strong_classical", "lean_classical", "lean_basal_like", "strong_basal_like"),
  right = TRUE
)

labels <- tibble(
  sample = tumor_samples,
  schema = "PurIST_Moffitt_tumor_intrinsic",
  purist_linear_score = linear_score,
  purist_probability_basal = prob_basal,
  purist_subtype = subtype,
  purist_graded_call = as.character(graded)
) |>
  bind_cols(as_tibble(pair_mat))

# This file records sample-level labels computed by the fixed published PurIST classifier.
write_tsv(labels, "resources/pdac_subtype_labels.tsv")
write_tsv(labels, file.path(out_dir, "pdac_subtype_labels_purist.tsv"))

if (!"sample" %in% names(scores)) {
  stop("TCGA module score table must contain sample.", call. = FALSE)
}
module_cols <- setdiff(names(scores), "sample")
df <- inner_join(labels, scores, by = "sample")
if (nrow(df) < 50L) {
  stop("Too few samples overlap PurIST calls and module scores: ", nrow(df), call. = FALSE)
}

benchmark <- bind_rows(lapply(module_cols, function(m) {
  x <- suppressWarnings(as.numeric(df[[m]]))
  basal <- df$purist_subtype == "basal_like"
  n_b <- sum(basal & is.finite(x))
  n_c <- sum(!basal & is.finite(x))
  p <- if (n_b >= 5L && n_c >= 5L) {
    tryCatch(
      stats::wilcox.test(x[basal], x[!basal], exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  } else NA_real_
  auc <- rank_auc(x, basal)
  rho <- suppressWarnings(
    stats::cor(
      x,
      df$purist_probability_basal,
      method = "spearman",
      use = "complete.obs"
    )
  )
  tibble(
    label = "PurIST_basal_like_vs_classical",
    module = m,
    n = sum(is.finite(x)),
    n_classes = 2L,
    n_basal_like = n_b,
    n_classical = n_c,
    basal_median = median(x[basal], na.rm = TRUE),
    classical_median = median(x[!basal], na.rm = TRUE),
    basal_minus_classical_median = basal_median - classical_median,
    auc_basal = auc,
    absolute_auc = ifelse(is.finite(auc), max(auc, 1 - auc), NA_real_),
    cliffs_delta_basal = ifelse(is.finite(auc), 2 * auc - 1, NA_real_),
    spearman_rho_with_basal_probability = rho,
    kruskal_p = p,
    mean_by_class = paste0(
      "basal_like=", signif(mean(x[basal], na.rm = TRUE), 4),
      ";classical=", signif(mean(x[!basal], na.rm = TRUE), 4)
    )
  )
})) |>
  mutate(fdr = p.adjust(.data$kruskal_p, method = "BH")) |>
  arrange(.data$fdr, desc(.data$absolute_auc))

write_tsv(benchmark, file.path(out_dir, "pdac_subtype_module_benchmark.tsv"))

summary_tbl <- tibble(
  method = "PurIST published 8-TSP single-sample classifier",
  expression_input_source = input_source,
  expression_input_file = expr_file,
  n_tcga_paad_tumors = nrow(labels),
  n_basal_like = sum(labels$purist_subtype == "basal_like"),
  n_classical = sum(labels$purist_subtype == "classical"),
  fraction_basal_like = mean(labels$purist_subtype == "basal_like"),
  n_modules_benchmarked = nrow(benchmark),
  min_fdr = suppressWarnings(min(benchmark$fdr, na.rm = TRUE)),
  max_absolute_auc = suppressWarnings(max(benchmark$absolute_auc, na.rm = TRUE))
)
write_tsv(summary_tbl, file.path(out_dir, "pdac_subtype_module_benchmark_summary.tsv"))

cat(
  "PurIST subtype benchmark complete: ", nrow(labels), " tumors; ",
  sum(labels$purist_subtype == "basal_like"), " basal-like; ",
  sum(labels$purist_subtype == "classical"), " classical; ",
  nrow(benchmark), " module tests.\n",
  sep = ""
)
