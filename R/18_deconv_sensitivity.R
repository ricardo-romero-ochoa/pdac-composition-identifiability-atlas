# Deconvolution sensitivity layer.
# This module compares the original atlas TME classification with independent
# marker-score deconvolution surrogates and optional package-based scores.

get_deconv_sensitivity_cfg <- function(cfg) get_validation_subcfg(cfg, "deconv_sensitivity")

builtin_deconv_marker_sets <- function() {
  list(
    caf_fibroblast = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL12A1", "LUM", "DCN", "FAP", "ACTA2", "TAGLN", "PDGFRB", "THY1", "MMP2"),
    immune_myeloid = c("PTPRC", "LST1", "TYROBP", "AIF1", "CD68", "CD14", "LYZ", "FCER1G", "ITGAM", "CSF1R"),
    t_cell = c("CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CCR7", "CD8A", "CD8B", "NKG7", "GZMB"),
    b_cell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "MZB1", "JCHAIN"),
    endothelial = c("PECAM1", "VWF", "CLDN5", "KDR", "ENG", "ESAM", "RAMP2"),
    ductal_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "SOX9", "KRT7"),
    acinar_normal = c("PRSS1", "PRSS2", "CPA1", "CPB1", "CTRB1", "CTRB2", "AMY2A", "REG1A", "CELA3A"),
    endocrine_normal = c("INS", "GCG", "SST", "PPY", "CHGA", "CHGB", "IAPP")
  )
}

compute_marker_deconv_scores <- function(expr, marker_sets = builtin_deconv_marker_sets(), min_genes = 3) {
  expr <- standardize_gene_matrix(expr)
  z <- zscore_rows_safe(expr)
  scores <- purrr::map(marker_sets, function(genes) {
    g <- intersect(normalise_gene_symbols(genes), rownames(z))
    if (length(g) < min_genes) return(setNames(rep(NA_real_, ncol(z)), colnames(z)))
    colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  })
  out <- as.data.frame(scores, check.names = FALSE)
  out$sample <- colnames(z)
  out |> dplyr::relocate(.data$sample)
}

compute_optional_gsva_deconv_scores <- function(expr, marker_sets = builtin_deconv_marker_sets()) {
  expr <- standardize_gene_matrix(expr)
  sets <- lapply(marker_sets, function(g) intersect(normalise_gene_symbols(g), rownames(expr)))
  sets <- sets[lengths(sets) >= 3]
  if (!length(sets) || !requireNamespace("GSVA", quietly = TRUE)) return(NULL)
  mat <- tryCatch({
    if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
      par <- GSVA::ssgseaParam(expr, sets, normalize = TRUE)
      GSVA::gsva(par)
    } else {
      GSVA::gsva(expr, sets, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(mat)) return(NULL)
  as.data.frame(t(mat), check.names = FALSE) |> tibble::rownames_to_column("sample")
}

summarise_deconv_condition_effects <- function(score_tbl, meta, dataset, group_col = "condition") {
  if (!"sample" %in% colnames(meta)) meta$sample <- rownames(meta)
  purrr::map_dfr(setdiff(colnames(score_tbl), "sample"), function(s) {
    vec <- score_tbl[[s]]; names(vec) <- score_tbl$sample
    score_test_table(vec, meta, group_col = group_col, positive = "tumor", dataset = dataset, module = s)
  })
}

gene_score_correlations <- function(expr, score_tbl, dataset, score_cols = NULL, max_genes = NULL) {
  expr <- standardize_gene_matrix(expr)
  common <- intersect(colnames(expr), score_tbl$sample)
  if (length(common) < 5) return(tibble::tibble())
  expr <- expr[, common, drop = FALSE]
  score_tbl <- score_tbl |> dplyr::filter(.data$sample %in% common) |> dplyr::arrange(match(.data$sample, common))
  score_cols <- score_cols %||% setdiff(colnames(score_tbl), "sample")
  if (!is.null(max_genes) && nrow(expr) > max_genes) {
    vars <- matrixStats::rowVars(expr, na.rm = TRUE)
    expr <- expr[order(vars, decreasing = TRUE)[seq_len(max_genes)], , drop = FALSE]
  }
  purrr::map_dfr(score_cols, function(sc) {
    v <- score_tbl[[sc]]
    if (sum(is.finite(v)) < 5 || stats::sd(v, na.rm = TRUE) == 0) return(tibble::tibble())
    r <- apply(expr, 1, function(g) suppressWarnings(stats::cor(g, v, method = "spearman", use = "pairwise.complete.obs")))
    tibble::tibble(dataset = dataset, gene = rownames(expr), score = sc, rho = as.numeric(r))
  })
}

meta_summarise_correlations <- function(cor_tbl) {
  if (!nrow(cor_tbl)) return(tibble::tibble())
  cor_tbl |>
    dplyr::group_by(.data$gene, .data$score) |>
    dplyr::summarise(k = dplyr::n(), median_rho = stats::median(.data$rho, na.rm = TRUE), direction_concordance = mean(sign(.data$rho) == sign(stats::median(.data$rho, na.rm = TRUE)), na.rm = TRUE), .groups = "drop")
}

run_deconv_sensitivity <- function(dataset_list, cfg) {
  dcfg <- get_deconv_sensitivity_cfg(cfg)
  marker_sets <- builtin_deconv_marker_sets()
  max_genes_cor <- dcfg$max_genes_for_correlation %||% 6000
  score_outputs <- purrr::imap(dataset_list, function(ds, dsid) {
    prepared <- prepare_dataset_for_downstream_scores(ds)
    marker <- compute_marker_deconv_scores(prepared$expr, marker_sets = marker_sets, min_genes = dcfg$min_marker_genes %||% 3) |>
      dplyr::mutate(method = "marker_zscore", .before = 1)
    gsva <- compute_optional_gsva_deconv_scores(prepared$expr, marker_sets = marker_sets)
    if (!is.null(gsva)) gsva <- gsva |> dplyr::mutate(method = "gsva_ssgsea", .before = 1)
    list(dataset = dsid, meta = prepared$meta, expr = prepared$expr, scores = dplyr::bind_rows(marker, gsva))
  })
  score_long <- purrr::map_dfr(score_outputs, function(x) x$scores |> dplyr::mutate(dataset = x$dataset, .before = 1))
  condition_effects <- purrr::map_dfr(score_outputs, function(x) {
    split(x$scores, x$scores$method) |>
      purrr::imap_dfr(function(tbl, method) {
        tbl2 <- tbl |> dplyr::select(-method)
        summarise_deconv_condition_effects(tbl2, x$meta, x$dataset) |> dplyr::mutate(method = method, .before = 2)
      })
  }) |>
    dplyr::mutate(fdr = p.adjust(.data$p, method = "BH"))
  cors <- purrr::map_dfr(score_outputs, function(x) {
    split(x$scores, x$scores$method) |>
      purrr::imap_dfr(function(tbl, method) {
        tbl2 <- tbl |> dplyr::select(-method)
        gene_score_correlations(x$expr, tbl2, x$dataset, max_genes = max_genes_cor) |> dplyr::mutate(method = method, .before = 2)
      })
  })
  cor_meta <- cors |> dplyr::group_by(.data$method, .data$gene, .data$score) |>
    dplyr::summarise(k = dplyr::n(), median_rho = stats::median(.data$rho, na.rm = TRUE), direction_concordance = mean(sign(.data$rho) == sign(stats::median(.data$rho, na.rm = TRUE)), na.rm = TRUE), .groups = "drop")
  tme_classes_file <- file.path(cfg$project$output_dir, "tables", "tumor_tme_classified_programs.tsv")
  if (file.exists(tme_classes_file)) {
    tme <- readr::read_tsv(tme_classes_file, show_col_types = FALSE)
    gcol <- intersect(c("feature", "gene"), colnames(tme))[1]
    tme <- tme |> dplyr::mutate(gene = normalise_gene_symbols(.data[[gcol]]))
    overlap <- cor_meta |> dplyr::left_join(tme |> dplyr::select("gene", dplyr::everything()), by = "gene") |>
      dplyr::group_by(.data$method, .data$score, .data$tme_class) |>
      dplyr::summarise(n = dplyr::n(), median_abs_rho = stats::median(abs(.data$median_rho), na.rm = TRUE), .groups = "drop")
  } else {
    overlap <- tibble::tibble()
  }
  write_tsv_safe(score_long, file.path(external_table_dir(cfg), "deconv_sensitivity_scores.tsv.gz"))
  write_tsv_safe(condition_effects, file.path(external_table_dir(cfg), "deconv_sensitivity_condition_effects.tsv"))
  write_tsv_safe(cor_meta, file.path(external_table_dir(cfg), "deconv_sensitivity_gene_correlations_meta.tsv.gz"))
  write_tsv_safe(overlap, file.path(external_table_dir(cfg), "deconv_sensitivity_tme_class_overlap.tsv"))
  if (nrow(condition_effects)) {
    p <- condition_effects |>
      dplyr::filter(is.finite(.data$logFC)) |>
      ggplot2::ggplot(ggplot2::aes(x = module, y = logFC)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(method ~ dataset, scales = "free_x") +
      ggplot2::labs(x = "Deconvolution signature", y = "Tumor-control score shift", title = "Deconvolution sensitivity: condition effects") +
      ggplot2::theme_bw(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
    ggplot2::ggsave(file.path(external_figure_dir(cfg), "deconv_sensitivity_condition_effects.png"), p, width = 13, height = 7, dpi = 300)
  }
  list(scores = score_long, condition_effects = condition_effects, correlation_meta = cor_meta, overlap = overlap)
}
