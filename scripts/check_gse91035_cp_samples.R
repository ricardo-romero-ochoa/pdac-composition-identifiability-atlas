#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

path <- "results/tables/curated_metadata_all.tsv"
if (!file.exists(path)) {
  stop("Missing ", path, ". Run: Rscript scripts/export_curated_metadata.R", call. = FALSE)
}

meta <- readr::read_tsv(path, show_col_types = FALSE)
cp <- meta |>
  filter(dataset == "GSE91035", sample %in% c("GSM2420007", "GSM2420010")) |>
  select(any_of(c("dataset", "sample", "title", "condition", "include", "notes")))

print(cp, n = Inf)

if (nrow(cp) != 2 || !all(cp$condition == "cp", na.rm = TRUE) || any(cp$include %in% TRUE, na.rm = TRUE)) {
  stop("CP sample check failed: GSM2420007/GSM2420010 must be condition='cp' and include=FALSE.", call. = FALSE)
}

message("GSE91035 CP sample check passed.")
