#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

path <- "results/tables/curated_metadata_all.tsv"
if (!file.exists(path)) {
  stop("Missing ", path, ". Run: Rscript scripts/export_curated_metadata.R", call. = FALSE)
}

meta <- readr::read_tsv(path, show_col_types = FALSE)
if (!all(c("dataset", "condition", "include") %in% names(meta))) {
  stop("curated_metadata_all.tsv is missing one of: dataset, condition, include", call. = FALSE)
}

counts <- meta |>
  filter(isTRUE(include) | include == TRUE) |>
  count(dataset, condition, name = "n") |>
  tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) |>
  arrange(dataset)

cat("Included samples used by the primary pipeline:\n")
print(counts, n = Inf)

# Confirmed manual GEO audit as of v2.4:
# - GSE71989 has 13 PDAC/tumor and 8 control samples.
# - GSE91035 has 25 PDAC/tumor, 8 control, 15 benign, and two CP samples
#   (GSM2420007, GSM2420010). CP is documented in metadata but excluded from
#   primary analyses by default.

all_counts <- meta |>
  count(dataset, condition, name = "n_all") |>
  tidyr::pivot_wider(names_from = condition, values_from = n_all, values_fill = 0) |>
  arrange(dataset)

cat("\nAll curated samples, including excluded CP/unlabeled samples if present:\n")
print(all_counts, n = Inf)

required_two_level <- c("GSE15471", "GSE28735", "GSE62165", "GSE16515", "GSE71989")
for (ds in required_two_level) {
  x <- meta |> filter(dataset == ds, include == TRUE)
  lev <- sort(unique(x$condition[!is.na(x$condition)]))
  if (!all(c("control", "tumor") %in% lev)) {
    stop(ds, " does not contain both control and tumor after curation. Observed: ", paste(lev, collapse = ", "), call. = FALSE)
  }
}

x <- meta |> filter(dataset == "GSE91035", include == TRUE)
lev <- sort(unique(x$condition[!is.na(x$condition)]))
if (!all(c("control", "benign", "tumor") %in% lev)) {
  stop("GSE91035 does not contain control, benign, and tumor after curation. Observed: ", paste(lev, collapse = ", "), call. = FALSE)
}
if ("cp" %in% lev) {
  stop("GSE91035 CP samples are included in the primary analysis. They should be documented but excluded by default.", call. = FALSE)
}
cp_samples <- meta |> filter(dataset == "GSE91035", sample %in% c("GSM2420007", "GSM2420010"))
if (nrow(cp_samples) != 2 || !all(cp_samples$condition == "cp", na.rm = TRUE)) {
  stop("GSE91035 CP samples GSM2420007/GSM2420010 were not curated as condition == 'cp'.", call. = FALSE)
}
if (any(cp_samples$include %in% TRUE, na.rm = TRUE)) {
  stop("GSE91035 CP samples GSM2420007/GSM2420010 should be include == FALSE by default.", call. = FALSE)
}

count_value <- function(tbl, label) {
  val <- tbl$n[tbl$condition == label]
  if (!length(val)) 0L else as.integer(val[[1]])
}

gse71989 <- meta |> filter(dataset == "GSE71989", include == TRUE) |> count(condition)
if (count_value(gse71989, "tumor") != 13L || count_value(gse71989, "control") != 8L) {
  stop("GSE71989 expected included counts are 13 tumor and 8 control after GEO audit.", call. = FALSE)
}

gse91035 <- meta |> filter(dataset == "GSE91035", include == TRUE) |> count(condition)
if (count_value(gse91035, "tumor") != 25L ||
    count_value(gse91035, "control") != 8L ||
    count_value(gse91035, "benign") != 15L) {
  stop("GSE91035 expected included counts are 25 tumor, 8 control, and 15 benign; CP excluded.", call. = FALSE)
}

message("Metadata condition checks passed, including v2.4 CP handling.")
