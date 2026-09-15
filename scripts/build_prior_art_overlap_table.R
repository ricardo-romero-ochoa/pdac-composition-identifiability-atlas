#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})
source("R/00_utils.R")

cfg <- read_config("config/atlas_config.yml")
out_dir <- file.path(cfg$project$output_dir %||% "results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prior_file <- "resources/prior_art_pdac_meta_analysis_gene_lists.tsv"
tier1_file <- file.path(out_dir, "program_tier1_core_signature.tsv")
if (!file.exists(tier1_file)) stop("Missing ", tier1_file, ". Run the atlas pipeline first.", call. = FALSE)
if (!file.exists(prior_file)) stop("Missing ", prior_file, ". Run scripts/create_prior_art_resource_from_literature.R or provide curated file.", call. = FALSE)

prior <- read_tsv(prior_file, col_types = cols(.default = col_character()), show_col_types = FALSE) |>
  mutate(gene = normalise_gene_symbols(.data$gene)) |>
  filter(!is.na(.data$gene), nzchar(.data$gene))
if (!nrow(prior) || any(prior$study_id == "example_prior_meta_analysis")) {
  stop("Prior-art file is empty or still contains example rows. Curate published PDAC meta-analysis gene lists before running.", call. = FALSE)
}

required_cols <- c("study_id", "list_name", "gene")
missing_cols <- setdiff(required_cols, names(prior))
if (length(missing_cols)) stop("Prior-art file missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

tier1 <- read_tsv(tier1_file, show_col_types = FALSE)
gene_col <- intersect(c("feature", "gene", "symbol"), names(tier1))[1]
if (is.na(gene_col)) stop("Tier 1 table has no feature/gene/symbol column.", call. = FALSE)
tier1_genes <- unique(normalise_gene_symbols(tier1[[gene_col]]))
tier1_genes <- tier1_genes[!is.na(tier1_genes) & nzchar(tier1_genes)]

by_list <- prior |>
  group_by(.data$study_id, .data$list_name) |>
  summarise(
    prior_genes = n_distinct(.data$gene),
    overlap_genes = n_distinct(.data$gene[.data$gene %in% tier1_genes]),
    fraction_prior_recovered = overlap_genes / prior_genes,
    fraction_tier1_explained = overlap_genes / length(tier1_genes),
    jaccard = overlap_genes / n_distinct(c(.data$gene, tier1_genes)),
    overlapping_gene_list = paste(sort(unique(.data$gene[.data$gene %in% tier1_genes])), collapse = ";"),
    .groups = "drop"
  ) |>
  arrange(desc(.data$fraction_prior_recovered), desc(.data$overlap_genes))

# Study-level metadata: keep the first non-missing value for each descriptive field.
# These are identifiers/descriptors, so force a character return type across all
# groups. This avoids dplyr/vctrs failures when an identifier such as PMID is
# present for some studies but missing for others.
desc_cols <- intersect(c("citation_key", "year", "title", "source_type", "datasets_used", "reusable_gene_lists", "extraction_status", "source_url", "pmid", "doi", "notes"), names(prior))

first_nonmissing_chr <- function(x) {
  y <- as.character(x)
  idx <- which(!is.na(y) & nzchar(trimws(y)))
  if (!length(idx)) return(NA_character_)
  y[[idx[[1]]]]
}

study_meta <- prior |>
  group_by(.data$study_id) |>
  summarise(across(all_of(desc_cols), first_nonmissing_chr), .groups = "drop")

summary <- prior |>
  group_by(.data$study_id) |>
  summarise(
    curated_gene_lists = n_distinct(.data$list_name),
    curated_genes = n_distinct(.data$gene),
    overlap_genes = n_distinct(.data$gene[.data$gene %in% tier1_genes]),
    fraction_prior_recovered = overlap_genes / curated_genes,
    fraction_tier1_explained = overlap_genes / length(tier1_genes),
    overlapping_gene_list = paste(sort(unique(.data$gene[.data$gene %in% tier1_genes])), collapse = ";"),
    .groups = "drop"
  ) |>
  left_join(study_meta, by = "study_id") |>
  arrange(desc(.data$fraction_prior_recovered), desc(.data$overlap_genes))

write_tsv(prior, file.path(out_dir, "prior_art_curated_gene_lists_input.tsv"))
write_tsv(by_list, file.path(out_dir, "prior_art_overlap_by_gene_list.tsv"))
write_tsv(summary, file.path(out_dir, "prior_art_overlap_summary.tsv"))
write_tsv(tibble(tier1_gene = tier1_genes), file.path(out_dir, "prior_art_overlap_tier1_universe.tsv"))
cat("Wrote prior-art overlap outputs for ", nrow(summary), " studies and ", nrow(by_list), " gene lists.\n", sep = "")
