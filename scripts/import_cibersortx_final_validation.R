#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})
source("R/00_utils.R")

cfg <- read_config("config/atlas_config.yml")
result_dir <- "results/external_deconvolution_inputs/cibersortx/results"
out_ext <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
ref_dir <- "data/external/reference_deconvolution"
dir.create(out_ext, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE15471", "GSE28735", "GSE62165", "GSE16515", "GSE71989", "GSE91035")
files <- file.path(result_dir, paste0(datasets, "_CIBERSORTx_Adjusted.txt"))
missing <- files[!file.exists(files) | is.na(file.info(files)$size) | file.info(files)$size == 0]
if (length(missing)) {
  stop(
    "Missing CIBERSORTx result files:\n  ",
    paste(missing, collapse = "\n  "),
    "\nRun scripts/run_cibersortx_final_validation.sh first.",
    call. = FALSE
  )
}

canonical_compartment <- function(x) {
  z <- tolower(gsub("[^a-z0-9]+", "_", x))
  case_when(
    grepl("malignan|cancer|tumou?r", z) ~ "malignant_epithelial",
    grepl("duct", z) ~ "ductal",
    grepl("acinar", z) ~ "acinar",
    grepl("fibro|caf|stellate|pericyte|vsmc|strom|mesench", z) ~ "fibroblast_stromal",
    grepl("immune|myeloid|macroph|monocyte|lymph|t_cell|b_cell|nk|mast|dendritic|neutroph|plasma", z) ~ "immune",
    grepl("endothel|vascular|capillary|arterial|venous", z) ~ "endothelial",
    grepl("endocrine|alpha|beta|delta|gamma|islet", z) ~ "endocrine",
    TRUE ~ NA_character_
  )
}

read_one <- function(dataset, path) {
  x <- read_tsv(path, show_col_types = FALSE)
  sample_col <- intersect(c("Mixture", "mixture", "Sample", "sample"), names(x))[1]
  if (is.na(sample_col)) {
    stop(dataset, ": cannot identify CIBERSORTx mixture/sample column in ", path, call. = FALSE)
  }

  qc_names <- intersect(
    c("P-value", "P.value", "Pvalue", "Correlation", "RMSE", "Absolute score"),
    names(x)
  )
  candidate <- setdiff(names(x), c(sample_col, qc_names))
  mapped <- canonical_compartment(candidate)
  class_cols <- candidate[!is.na(mapped)]
  class_map <- tibble(
    source_column = class_cols,
    compartment = canonical_compartment(class_cols)
  ) |>
    distinct(.data$source_column, .keep_all = TRUE)

  required_classes <- c(
    "malignant_epithelial", "ductal", "acinar",
    "fibroblast_stromal", "immune", "endothelial", "endocrine"
  )
  absent <- setdiff(required_classes, unique(class_map$compartment))
  if (length(absent)) {
    stop(
      dataset, ": returned CIBERSORTx table does not contain all 7 expected classes. ",
      "Missing: ", paste(absent, collapse = ", "),
      ". Columns: ", paste(names(x), collapse = ", "),
      call. = FALSE
    )
  }

  prop <- x |>
    transmute(sample = as.character(.data[[sample_col]]), across(all_of(class_cols))) |>
    pivot_longer(-"sample", names_to = "source_column", values_to = "proportion") |>
    left_join(class_map, by = "source_column") |>
    mutate(
      dataset = dataset,
      method = "CIBERSORTx_Smode_PDac7_scRNA",
      reference = "PDAC_scRNA_7class",
      proportion = suppressWarnings(as.numeric(.data$proportion)),
      .before = 1
    ) |>
    filter(is.finite(.data$proportion)) |>
    select(.data$method, .data$reference, .data$dataset, .data$sample,
           .data$compartment, .data$proportion, .data$source_column)

  qc <- tibble(dataset = dataset, sample = as.character(x[[sample_col]]))
  for (nm in qc_names) {
    qc[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  }
  list(prop = prop, qc = qc, class_map = class_map)
}

res <- Map(read_one, datasets, files)
prop <- bind_rows(lapply(res, `[[`, "prop"))
qc <- bind_rows(lapply(res, `[[`, "qc"))
class_map <- bind_rows(Map(
  function(ds, obj) mutate(obj$class_map, dataset = ds, .before = 1),
  datasets, res
))

write_tsv(prop, file.path(ref_dir, "cibersortx_proportions.tsv"))
write_tsv(qc, file.path(out_ext, "reference_deconvolution_qc_samples.tsv"))
write_tsv(class_map, file.path(out_ext, "cibersortx_returned_class_map.tsv"))

# Make CIBERSORTx the canonical final reference-based result. Preserve the old
# local NNLS file separately if it exists.
canonical_file <- file.path(ref_dir, "reference_deconvolution_proportions.tsv")
if (file.exists(canonical_file)) {
  old <- read_tsv(canonical_file, show_col_types = FALSE)
  old_method <- if ("method" %in% names(old)) unique(as.character(old$method)) else character()
  if (!any(grepl("CIBERSORTx", old_method, fixed = TRUE))) {
    backup <- file.path(ref_dir, "reference_deconvolution_proportions_local_rankNNLS.tsv")
    if (!file.exists(backup)) file.copy(canonical_file, backup, overwrite = FALSE)
  }
}
write_tsv(prop, canonical_file)

marker_file <- "results/tables/deconvolution_scores_long.tsv"
if (!file.exists(marker_file)) {
  stop("Missing atlas marker-score table: ", marker_file, call. = FALSE)
}
marker <- read_tsv(marker_file, show_col_types = FALSE)
if (!all(c("dataset", "sample", "score", "value") %in% names(marker))) {
  stop(
    "Expected marker table schema dataset/sample/score/value in ", marker_file,
    ". Found: ", paste(names(marker), collapse = ", "),
    call. = FALSE
  )
}
marker2 <- marker |>
  transmute(
    dataset = as.character(.data$dataset),
    sample = as.character(.data$sample),
    marker_signature = as.character(.data$score),
    marker_score = suppressWarnings(as.numeric(.data$value))
  ) |>
  filter(is.finite(.data$marker_score))

joined <- inner_join(
  prop |> select(-.data$source_column),
  marker2,
  by = c("dataset", "sample")
)

expected_pair_flag <- function(compartment, marker_signature) {
  (compartment == "fibroblast_stromal" & marker_signature == "stromal_caf") |
  (compartment == "immune" & marker_signature == "immune_pan") |
  (compartment == "endothelial" & marker_signature == "endothelial") |
  (compartment == "acinar" & marker_signature == "acinar_pancreas") |
  (compartment == "ductal" & marker_signature == "ductal_epithelial")
}

safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 5 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  z <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(z$estimate), p = z$p.value)
}

summary <- joined |>
  group_by(.data$dataset, .data$method, .data$compartment, .data$marker_signature) |>
  group_modify(~ {
    z <- safe_cor_test(.x$proportion, .x$marker_score)
    tibble(
      n = sum(is.finite(.x$proportion) & is.finite(.x$marker_score)),
      spearman_rho = z[["rho"]],
      spearman_p = z[["p"]]
    )
  }) |>
  ungroup() |>
  mutate(
    expected_pair = expected_pair_flag(.data$compartment, .data$marker_signature)
  ) |>
  arrange(desc(.data$expected_pair), .data$compartment, .data$dataset)

write_tsv(summary, file.path(out_ext, "reference_deconvolution_summary.tsv"))
write_tsv(joined, file.path(out_ext, "reference_deconvolution_marker_score_joined.tsv"))

expected <- summary |>
  filter(.data$expected_pair) |>
  mutate(fdr_within_expected = p.adjust(.data$spearman_p, method = "BH"))
write_tsv(expected, file.path(out_ext, "reference_deconvolution_expected_pair_by_dataset.tsv"))

clamp_rho <- function(x) pmin(pmax(x, -0.999999), 0.999999)

meta <- expected |>
  group_by(.data$compartment, .data$marker_signature) |>
  summarise(
    n_datasets = sum(is.finite(.data$spearman_rho)),
    median_rho = median(.data$spearman_rho, na.rm = TRUE),
    min_rho = min(.data$spearman_rho, na.rm = TRUE),
    max_rho = max(.data$spearman_rho, na.rm = TRUE),
    positive_fraction = mean(.data$spearman_rho > 0, na.rm = TRUE),
    fisher_z_weighted_rho = {
      ok <- is.finite(.data$spearman_rho) & .data$n > 3
      if (!any(ok)) NA_real_ else {
        w <- pmax(.data$n[ok] - 3, 1)
        tanh(sum(w * atanh(clamp_rho(.data$spearman_rho[ok]))) / sum(w))
      }
    },
    .groups = "drop"
  ) |>
  arrange(desc(.data$median_rho))

write_tsv(meta, file.path(out_ext, "reference_deconvolution_expected_pair_meta.tsv"))

# QC aggregation from CIBERSORTx returned columns.
qc_summary <- qc |>
  group_by(.data$dataset) |>
  summarise(
    n_samples = n(),
    median_cibersortx_correlation = if ("Correlation" %in% names(qc)) median(.data$Correlation, na.rm = TRUE) else NA_real_,
    median_cibersortx_rmse = if ("RMSE" %in% names(qc)) median(.data$RMSE, na.rm = TRUE) else NA_real_,
    fraction_p_le_005 = if ("P-value" %in% names(qc)) mean(.data$`P-value` <= 0.05, na.rm = TRUE) else
      if ("P.value" %in% names(qc)) mean(.data$P.value <= 0.05, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )
write_tsv(qc_summary, file.path(out_ext, "reference_deconvolution_qc.tsv"))

cat("\nCIBERSORTx fractions imported: ", nrow(prop), " sample-class rows.\n", sep = "")
cat("Expected-pair within-cohort benchmark:\n")
print(expected |>
  select(.data$dataset, .data$compartment, .data$marker_signature,
         .data$n, .data$spearman_rho, .data$fdr_within_expected))
cat("\nAcross-cohort expected-pair summary:\n")
print(meta)
cat("\nWrote canonical reference-based summary:\n  ",
    file.path(out_ext, "reference_deconvolution_summary.tsv"), "\n", sep = "")
