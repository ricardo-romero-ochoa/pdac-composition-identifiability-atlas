# Microarray preprocessing and gene mapping ----------------------------------

symbol_from_featuredata <- function(es) {
  fd <- Biobase::fData(es) |> as.data.frame()
  candidates <- c(
    "Gene Symbol", "GENE_SYMBOL", "Gene symbol", "gene_symbol", "Symbol", "SYMBOL",
    "Gene.symbol", "GENE_SYMBOL_ASSIGNMENT", "gene_assignment", "Gene Assignment",
    "ILMN_Gene", "GB_ACC", "ORF", "gene" 
  )
  col <- first_present_col(fd, candidates)
  if (!is.na(col)) return(clean_symbol(fd[[col]]))

  # Search flexibly for symbol-like columns.
  symbol_cols <- grep("symbol|gene.*assignment|hgnc", names(fd), ignore.case = TRUE, value = TRUE)
  for (cc in symbol_cols) {
    s <- clean_symbol(fd[[cc]])
    if (sum(!is.na(s)) > 1000) return(s)
  }
  rep(NA_character_, nrow(fd))
}

symbol_from_annotation_package <- function(es) {
  platform <- Biobase::annotation(es) %||% ""
  probes <- rownames(Biobase::exprs(es))
  if (platform == "GPL570" && quiet_require("hgu133plus2.db")) {
    return(AnnotationDbi::mapIds(hgu133plus2.db::hgu133plus2.db, keys = probes, keytype = "PROBEID", column = "SYMBOL", multiVals = "first") |> clean_symbol())
  }
  if (platform %in% c("GPL6244", "GPL5175") && quiet_require("hugene10sttranscriptcluster.db")) {
    return(AnnotationDbi::mapIds(hugene10sttranscriptcluster.db::hugene10sttranscriptcluster.db, keys = probes, keytype = "PROBEID", column = "SYMBOL", multiVals = "first") |> clean_symbol())
  }
  rep(NA_character_, length(probes))
}

preprocess_one <- function(gse_id, es, meta, ds_cfg, cfg) {
  cache <- file.path(cfg$project$data_dir, "processed", paste0(gse_id, "_gene_expression.rds"))
  if (file.exists(cache)) return(readRDS(cache))

  expr <- Biobase::exprs(es) |> as_numeric_matrix()
  if (isTRUE(cfg$analysis$log2_if_needed)) expr <- log2_if_needed(expr)
  if (isTRUE(cfg$analysis$normalize_between_arrays)) expr <- limma::normalizeBetweenArrays(expr, method = "quantile")

  symbols <- symbol_from_featuredata(es)
  if (sum(!is.na(symbols)) < cfg$analysis$min_genes_per_dataset) {
    symbols2 <- symbol_from_annotation_package(es)
    if (sum(!is.na(symbols2)) > sum(!is.na(symbols))) symbols <- symbols2
  }

  gene_expr <- collapse_by_symbol(expr, symbols, strategy = cfg$analysis$collapse_probe_strategy %||% "mean")
  keep_samples <- intersect(colnames(gene_expr), meta$sample[isTRUE(meta$include) | meta$include == TRUE])
  gene_expr <- gene_expr[, keep_samples, drop = FALSE]
  meta <- meta |> dplyr::filter(.data$sample %in% keep_samples) |> dplyr::arrange(match(.data$sample, keep_samples))
  gene_expr <- gene_expr[, meta$sample, drop = FALSE]

  obj <- list(
    dataset = gse_id,
    platform = Biobase::annotation(es),
    expr = gene_expr,
    meta = meta,
    config = ds_cfg
  )

  if (nrow(gene_expr) < cfg$analysis$min_genes_per_dataset) {
    warning(gse_id, ": only ", nrow(gene_expr), " genes mapped. Check GPL annotation; lncRNA-focused arrays may map fewer coding genes.")
  }
  saveRDS(obj, cache)
  obj
}

preprocess_all <- function(esets, metadata_list, cfg) {
  purrr::imap(esets, ~ preprocess_one(.y, .x, metadata_list[[.y]], cfg$datasets[[.y]], cfg))
}

common_genes <- function(dataset_list, min_datasets = NULL) {
  gs <- lapply(dataset_list, function(x) rownames(x$expr))
  if (is.null(min_datasets)) return(Reduce(intersect, gs))
  tab <- sort(table(unlist(gs)), decreasing = TRUE)
  names(tab)[tab >= min_datasets]
}

write_expression_summaries <- function(dataset_list, cfg) {
  summary <- purrr::map_dfr(dataset_list, function(ds) {
    ds$meta |>
      dplyr::count(.data$dataset, .data$condition, name = "n_samples") |>
      dplyr::mutate(platform = ds$platform, n_genes = nrow(ds$expr))
  })
  write_tsv(summary, file.path(cfg$project$output_dir, "tables", "dataset_sample_summary.tsv"))
  summary
}
