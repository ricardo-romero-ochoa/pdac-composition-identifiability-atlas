#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})
source("R/00_utils.R")

# This script validates the repository-bundled prior-art curation table. It does
# not scrape papers at run time; the curation is intentionally static and
# auditable. Update resources/prior_art_pdac_meta_analysis_gene_lists.tsv when
# adding new literature.
prior_file <- "resources/prior_art_pdac_meta_analysis_gene_lists.tsv"
if (!file.exists(prior_file)) stop("Missing ", prior_file, call. = FALSE)
prior <- read_tsv(prior_file, show_col_types = FALSE)
required <- c("study_id", "citation_key", "year", "title", "list_name", "gene", "source_url", "doi")
missing <- setdiff(required, names(prior))
if (length(missing)) stop("Prior-art curation missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
prior <- prior |> mutate(gene = normalise_gene_symbols(.data$gene)) |> filter(!is.na(.data$gene), nzchar(.data$gene))
if (!nrow(prior)) stop("Prior-art curation has no valid genes.", call. = FALSE)
if (any(prior$study_id == "example_prior_meta_analysis")) stop("Example prior-art row still present.", call. = FALSE)
summary <- prior |>
  count(.data$study_id, .data$citation_key, .data$year, .data$title, name = "n_curated_gene_rows") |>
  arrange(.data$year, .data$study_id)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write_tsv(prior, "results/tables/prior_art_curated_gene_lists_input.tsv")
write_tsv(summary, "results/tables/prior_art_curation_source_summary.tsv")
cat("Validated prior-art curation: ", n_distinct(prior$study_id), " studies, ", nrow(prior), " gene rows.\n", sep = "")
