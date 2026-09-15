#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})
source("R/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
strict <- !any(args %in% c("--allow-missing", "--no-fail"))
cfg <- read_config("config/atlas_config.yml")
out_dir <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_candidates <- c(
  "data/external/reference_deconvolution/cibersortx_proportions.tsv",
  "data/external/reference_deconvolution/bayesprism_proportions.tsv",
  "data/external/reference_deconvolution/reference_deconvolution_proportions.tsv"
)
input_file <- input_candidates[file.exists(input_candidates)][1]
marker_file <- "results/tables/deconvolution_scores_long.tsv"

if (is.na(input_file) || !file.exists(input_file)) {
  dir.create("data/external/reference_deconvolution", recursive = TRUE, showWarnings = FALSE)
  tmpl <- tibble(sample = character(), method = character(), compartment = character(), proportion = numeric(), reference = character())
  write_tsv(tmpl, "data/external/reference_deconvolution/reference_deconvolution_proportions_TEMPLATE.tsv")
  msg <- "Reference-based deconvolution proportions are missing. Provide CIBERSORTx/BayesPrism output at data/external/reference_deconvolution/reference_deconvolution_proportions.tsv."
  write_tsv(tibble(status = "missing", message = msg), file.path(out_dir, "reference_deconvolution_summary.tsv"))
  if (strict) stop(msg, call. = FALSE) else { cat(msg, "\n"); quit(status = 0) }
}
if (!file.exists(marker_file)) stop("Missing marker-average deconvolution scores: ", marker_file, call. = FALSE)

ref <- read_tsv(input_file, show_col_types = FALSE)
req <- c("sample", "compartment", "proportion")
miss <- setdiff(req, names(ref))
if (length(miss)) stop("Reference deconvolution file missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
if (!"method" %in% names(ref)) ref$method <- "reference_deconvolution"
ref <- ref |> mutate(sample = as.character(.data$sample), proportion = suppressWarnings(as.numeric(.data$proportion))) |> filter(is.finite(.data$proportion))

marker <- read_tsv(marker_file, show_col_types = FALSE)
# Expected marker table columns are dataset/sample/signature/score in most releases.
sample_col <- intersect(c("sample", "sample_id", "array", "gsm"), names(marker))[1]
# Native atlas schema: dataset / sample / score / value.
if (all(c("score", "value") %in% names(marker))) {
  sig_col <- "score"
  score_col <- "value"
} else {
  sig_col <- intersect(c("signature", "score_name", "cell_type", "method"), names(marker))[1]
  score_col <- intersect(c("value", "ssgsea_score", "zscore", "score"), names(marker))[1]
}
if (any(is.na(c(sample_col, sig_col, score_col)))) {
  stop("Cannot infer sample/signature/score columns in ", marker_file, call. = FALSE)
}
marker2 <- marker |>
  transmute(
    sample = as.character(.data[[sample_col]]),
    marker_signature = as.character(.data[[sig_col]]),
    marker_score = suppressWarnings(as.numeric(.data[[score_col]]))
  ) |>
  filter(is.finite(.data$marker_score))

joined <- inner_join(ref, marker2, by = "sample")
summary <- joined |>
  group_by(.data$method, .data$compartment, .data$marker_signature) |>
  summarise(
    n = n(),
    spearman_rho = suppressWarnings(
      cor(.data$proportion, .data$marker_score, method = "spearman", use = "complete.obs")
    ),
    .groups = "drop"
  ) |>
  mutate(
    expected_pair = case_when(
      .data$compartment == "fibroblast_stromal" & .data$marker_signature == "stromal_caf" ~ TRUE,
      .data$compartment == "immune" & .data$marker_signature == "immune_pan" ~ TRUE,
      .data$compartment == "endothelial" & .data$marker_signature == "endothelial" ~ TRUE,
      .data$compartment == "nonmalignant_epithelial" & .data$marker_signature %in% c("ductal_epithelial", "acinar_pancreas") ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  arrange(desc(.data$expected_pair), .data$method, .data$compartment, desc(abs(.data$spearman_rho)))

write_tsv(summary, file.path(out_dir, "reference_deconvolution_summary.tsv"))
write_tsv(joined, file.path(out_dir, "reference_deconvolution_marker_score_joined.tsv"))
cat("Wrote reference deconvolution benchmark from ", input_file, ": ", nrow(summary), " rows.\n", sep = "")
