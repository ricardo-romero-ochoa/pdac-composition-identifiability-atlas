#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
})
source("R/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
strict <- !any(args %in% c("--allow-partial", "--no-fail"))

cfg <- read_config("config/atlas_config.yml")
out_dir <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manual_file <- "resources/pdac_subtype_signature_gene_sets.tsv"
source_inventory <- "resources/pdac_subtype_signature_source_inventory.tsv"
cache_dir <- "data/external/subtype_signatures"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

module_files <- c(
  tier1_core = "results/tables/program_tier1_core_signature.tsv",
  tumor_associated_retained = "results/tables/program_tumor_intrinsic_retained_signature.tsv",
  composition_covarying = "results/tables/program_composition_covarying_signature.tsv"
)
read_module_genes <- function(path) {
  if (!file.exists(path)) return(character())
  x <- readr::read_tsv(path, show_col_types = FALSE)
  gene_col <- intersect(c("feature", "gene", "symbol"), names(x))[1]
  if (is.na(gene_col)) return(character())
  unique(normalise_gene_symbols(x[[gene_col]])) |> na.omit() |> as.character()
}
modules <- lapply(module_files, read_module_genes)
modules <- modules[lengths(modules) > 0]
if (!length(modules)) stop("No atlas module gene files found. Run the atlas pipeline first.", call. = FALSE)
universe <- unique(unlist(modules, use.names = FALSE))

# Helper for parsing small RData signature files from SignatureHeatmap. This is
# deliberately generic because those files are simple but can use different object names.
extract_gene_sets_from_object <- function(obj, source_id) {
  res <- list()
  add_vec <- function(v, nm) {
    vv <- unique(normalise_gene_symbols(as.character(v)))
    vv <- vv[!is.na(vv) & nzchar(vv)]
    if (length(vv) >= 3) res[[paste(source_id, nm, sep = "__")]] <<- vv
  }
  walk <- function(x, nm = "signature") {
    if (is.vector(x) && !is.list(x)) {
      add_vec(x, nm)
    } else if (is.data.frame(x)) {
      nms <- names(x)
      gene_col <- intersect(c("gene", "genes", "symbol", "Gene", "SYMBOL", "feature"), nms)[1]
      class_col <- intersect(c("subtype", "class", "signature", "group", "Subtype", "Class"), nms)[1]
      if (!is.na(gene_col) && !is.na(class_col)) {
        for (cl in unique(x[[class_col]])) add_vec(x[[gene_col]][x[[class_col]] == cl], paste(nm, cl, sep = "_"))
      } else if (!is.na(gene_col)) {
        add_vec(x[[gene_col]], nm)
      } else {
        for (cc in nms) if (is.character(x[[cc]])) add_vec(x[[cc]], paste(nm, cc, sep = "_"))
      }
    } else if (is.list(x)) {
      nms <- names(x); if (is.null(nms)) nms <- paste0("set", seq_along(x))
      for (i in seq_along(x)) walk(x[[i]], paste(nm, nms[[i]], sep = "_"))
    }
  }
  walk(obj, source_id)
  res
}

sig_tbl <- tibble()
# Manual/curated TSV first, if populated.
if (file.exists(manual_file)) {
  manual <- readr::read_tsv(manual_file, show_col_types = FALSE)
  if (nrow(manual) && all(c("signature_source", "subtype_class", "gene") %in% names(manual))) {
    sig_tbl <- manual |>
      transmute(signature_source = .data$signature_source, subtype_class = .data$subtype_class, gene = normalise_gene_symbols(.data$gene), direction = .data$direction %||% NA_character_) |>
      filter(!is.na(.data$gene), nzchar(.data$gene))
  }
}

# Fetch SignatureHeatmap RData sources listed in the source inventory. This gives
# a reproducible basal/classical benchmark when internet is available.
fetched_status <- tibble(source_id = character(), status = character(), n_sets = integer(), note = character())
if (file.exists(source_inventory)) {
  inv <- readr::read_tsv(source_inventory, show_col_types = FALSE)
  rdata_inv <- inv |> filter(grepl("RData", .data$source_type), !is.na(.data$url), nzchar(.data$url))
  for (i in seq_len(nrow(rdata_inv))) {
    sid <- rdata_inv$source_id[i]
    url <- rdata_inv$url[i]
    dest <- file.path(cache_dir, paste0(sid, ".RData"))
    status <- "not_attempted"; note <- ""; n_sets <- 0L
    try({
      if (!file.exists(dest)) utils::download.file(url, dest, mode = "wb", quiet = TRUE)
      env <- new.env(parent = emptyenv())
      loaded <- load(dest, envir = env)
      all_sets <- list()
      for (nm in loaded) all_sets <- c(all_sets, extract_gene_sets_from_object(env[[nm]], sid))
      if (length(all_sets)) {
        tmp <- bind_rows(lapply(names(all_sets), function(nm) {
          # Use the last token as subtype class where possible.
          subtype <- sub("^.*__(.*)$", "\\1", nm)
          tibble(signature_source = paste0("SignatureHeatmap_", sid), subtype_class = subtype, gene = all_sets[[nm]], direction = NA_character_)
        }))
        sig_tbl <- bind_rows(sig_tbl, tmp)
        n_sets <- length(all_sets)
        status <- "ok"
      } else {
        status <- "no_gene_sets_extracted"
      }
    }, silent = TRUE)
    fetched_status <- bind_rows(fetched_status, tibble(source_id = sid, status = status, n_sets = n_sets, note = note))
  }
}

sig_tbl <- sig_tbl |>
  mutate(gene = normalise_gene_symbols(.data$gene)) |>
  filter(!is.na(.data$gene), nzchar(.data$gene)) |>
  distinct(.data$signature_source, .data$subtype_class, .data$gene, .keep_all = TRUE)

write_tsv(fetched_status, file.path(out_dir, "pdac_subtype_signature_fetch_status.tsv"))
if (!nrow(sig_tbl)) {
  msg <- "No subtype signature gene sets available. Populate resources/pdac_subtype_signature_gene_sets.tsv or allow internet access for SignatureHeatmap downloads."
  write_tsv(tibble(status = "missing", message = msg), file.path(out_dir, "pdac_subtype_signature_overlap.tsv"))
  if (strict) stop(msg, call. = FALSE) else { cat(msg, "\n"); quit(status = 0) }
}
write_tsv(sig_tbl, file.path(out_dir, "pdac_subtype_signature_gene_sets_used.tsv"))

fisher_one <- function(a, b, bg) {
  ov <- length(intersect(a, b))
  mat <- matrix(c(ov, length(a) - ov, length(b) - ov, length(bg) - length(union(a, b))), nrow = 2)
  ft <- tryCatch(stats::fisher.test(mat, alternative = "greater"), error = function(e) NULL)
  list(or = if (is.null(ft)) NA_real_ else unname(ft$estimate), p = if (is.null(ft)) NA_real_ else ft$p.value)
}

res <- bind_rows(lapply(names(modules), function(mod) {
  a <- intersect(modules[[mod]], universe)
  sig_tbl |>
    group_by(.data$signature_source, .data$subtype_class) |>
    summarise(subtype_genes = list(unique(.data$gene)), .groups = "drop") |>
    rowwise() |>
    mutate(
      atlas_module = mod,
      atlas_genes = length(a),
      subtype_n_genes = length(.data$subtype_genes),
      overlap_genes = length(intersect(a, .data$subtype_genes)),
      fraction_subtype_recovered = overlap_genes / subtype_n_genes,
      fraction_atlas_explained = overlap_genes / atlas_genes,
      overlapping_gene_list = paste(sort(intersect(a, .data$subtype_genes)), collapse = ";"),
      fisher_or = fisher_one(a, .data$subtype_genes, universe)$or,
      fisher_p = fisher_one(a, .data$subtype_genes, universe)$p
    ) |>
    ungroup() |>
    select(.data$atlas_module, .data$signature_source, .data$subtype_class, .data$atlas_genes, .data$subtype_n_genes, .data$overlap_genes, .data$fraction_subtype_recovered, .data$fraction_atlas_explained, .data$fisher_or, .data$fisher_p, .data$overlapping_gene_list)
})) |>
  mutate(fdr = p.adjust(.data$fisher_p, method = "BH")) |>
  arrange(.data$fdr, desc(.data$overlap_genes))

write_tsv(res, file.path(out_dir, "pdac_subtype_signature_overlap.tsv"))
cat("Wrote subtype signature-overlap benchmark: ", nrow(res), " rows.\n", sep = "")
