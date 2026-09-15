# Publication-grade figures ---------------------------------------------------

save_plot <- function(p, filename, width = 7, height = 5, cfg = read_config()) {
  path <- file.path(cfg$project$output_dir, "figures", filename)
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 320, bg = "white")
  path
}

save_status_figure <- function(filename, message, cfg, width = 7, height = 4.5) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4, hjust = 0.5, vjust = 0.5) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = "Figure not generated")
  save_plot(p, filename, width = width, height = height, cfg = cfg)
}

safe_figure <- function(expr, filename, cfg, message = NULL) {
  tryCatch(
    expr,
    error = function(e) {
      msg <- message %||% paste("Skipped because:", conditionMessage(e))
      warning(msg, call. = FALSE)
      save_status_figure(filename, msg, cfg)
    }
  )
}

prepare_matrix_for_pca <- function(expr, min_fraction_finite = 0.8) {
  if (!nrow(expr) || !ncol(expr)) stop("Empty expression matrix supplied to PCA.")
  min_finite <- max(2, ceiling(ncol(expr) * min_fraction_finite))
  mat <- sanitize_numeric_matrix(
    expr,
    min_finite_per_row = min_finite,
    min_finite_per_col = 2,
    impute = "row_median",
    drop_zero_variance_rows = TRUE,
    drop_zero_variance_cols = FALSE
  )
  if (nrow(mat) < 2 || ncol(mat) < 3) {
    stop("Too few finite, variable genes/samples for PCA after filtering.")
  }
  t(mat)
}

prepare_matrix_for_heatmap <- function(mat, min_finite_per_row = 2, min_finite_per_col = 2) {
  mat <- sanitize_numeric_matrix(
    mat,
    min_finite_per_row = min_finite_per_row,
    min_finite_per_col = min_finite_per_col,
    impute = "row_median",
    drop_zero_variance_rows = TRUE,
    drop_zero_variance_cols = FALSE
  )
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Too few finite rows/columns for heatmap clustering.")
  }
  mat
}

plot_qc_pca_one <- function(ds, cfg) {
  safe_figure({
    prepared <- prepare_dataset_for_downstream_scores(ds)
    ds$expr <- prepared$expr
    ds$meta <- prepared$meta
    mat <- prepare_matrix_for_pca(ds$expr)
    pc <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
    pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
    df <- as.data.frame(pc$x[, 1:2, drop = FALSE]) |>
      tibble::rownames_to_column("sample") |>
      dplyr::left_join(ds$meta, by = "sample")
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$PC1, .data$PC2, color = .data$condition)) +
      ggplot2::geom_point(size = 2.2, alpha = 0.9) +
      ggplot2::labs(title = paste0(ds$dataset, " PCA QC"), x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"), color = "Group") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    save_plot(p, paste0("QC_PCA_", ds$dataset, ".png"), width = 5.5, height = 4.5, cfg = cfg)
  }, filename = paste0("QC_PCA_", ds$dataset, ".png"), cfg = cfg)
}

plot_all_qc <- function(dataset_list, cfg) {
  purrr::map_chr(dataset_list, plot_qc_pca_one, cfg = cfg)
}

plot_replication_heatmap <- function(meta_res, cfg, top_n = 100) {
  safe_figure({
    all <- meta_res$per_study
    top <- meta_res$meta |>
      dplyr::filter(is.finite(.data$meta_logFC)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::pull(.data$feature)
    if (!length(top)) stop("No finite top meta-analysis genes available for replication heatmap.")
    mat <- all |>
      dplyr::filter(.data$gene %in% top) |>
      dplyr::select("dataset", "gene", "logFC") |>
      dplyr::mutate(logFC = ifelse(is.finite(.data$logFC), .data$logFC, NA_real_)) |>
      tidyr::pivot_wider(names_from = "dataset", values_from = "logFC") |>
      tibble::column_to_rownames("gene") |>
      as.matrix()
    mat <- prepare_matrix_for_heatmap(mat, min_finite_per_row = 2, min_finite_per_col = 2)
    path <- file.path(cfg$project$output_dir, "figures", "replication_heatmap_top_meta_genes.png")
    grDevices::png(path, width = 2200, height = 2600, res = 300)
    pheatmap::pheatmap(
      mat,
      cluster_cols = ncol(mat) > 1,
      cluster_rows = nrow(mat) > 1,
      show_rownames = TRUE,
      fontsize_row = 5,
      main = "Per-study logFC for top meta-analysis genes"
    )
    grDevices::dev.off()
    path
  }, filename = "replication_heatmap_top_meta_genes.png", cfg = cfg)
}

plot_heterogeneity <- function(meta_res, cfg) {
  safe_figure({
    df <- meta_res$meta |>
      dplyr::filter(is.finite(.data$I2), is.finite(.data$meta_logFC)) |>
      dplyr::mutate(significant = .data$fdr < 0.05)
    if (!nrow(df)) stop("No finite meta-analysis rows available for heterogeneity plot.")
    p <- ggplot2::ggplot(df, ggplot2::aes(abs(.data$meta_logFC), .data$I2, alpha = .data$significant)) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::geom_hline(yintercept = 60, linetype = 2) +
      ggplot2::labs(x = "Absolute random-effects logFC", y = "I² heterogeneity (%)", title = "Effect size versus cross-study heterogeneity") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "none")
    save_plot(p, "heterogeneity_effect_size_scatter.png", width = 6, height = 4.5, cfg = cfg)
  }, filename = "heterogeneity_effect_size_scatter.png", cfg = cfg)
}

plot_forest_gene <- function(gene, meta_res, cfg) {
  safe_figure({
    df <- meta_res$per_study |>
      dplyr::filter(.data$gene == !!gene, is.finite(.data$logFC), is.finite(.data$se))
    if (!nrow(df)) stop(paste0("No finite per-study rows for gene ", gene, "."))
    pooled <- meta_res$meta |>
      dplyr::filter(.data$feature == !!gene, is.finite(.data$meta_logFC), is.finite(.data$ci_lb), is.finite(.data$ci_ub))
    df <- df |> dplyr::mutate(label = paste0(.data$dataset, " (n=", .data$n_tumor + .data$n_control, ")"))
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$logFC, stats::reorder(.data$label, .data$logFC))) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$logFC - 1.96 * .data$se, xmax = .data$logFC + 1.96 * .data$se), orientation = "y", height = 0.12) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "logFC", y = NULL, title = paste0("Forest plot: ", gene)) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    if (nrow(pooled)) {
      p <- p +
        ggplot2::geom_point(data = pooled, ggplot2::aes(x = .data$meta_logFC, y = "Pooled random effects"), inherit.aes = FALSE, size = 3, shape = 18) +
        ggplot2::geom_errorbar(data = pooled, ggplot2::aes(xmin = .data$ci_lb, xmax = .data$ci_ub, y = "Pooled random effects"), inherit.aes = FALSE, orientation = "y", height = 0.12)
    }
    save_plot(p, paste0("forest_gene_", make.names(gene), ".png"), width = 6.5, height = 4.5, cfg = cfg)
  }, filename = paste0("forest_gene_", make.names(gene), ".png"), cfg = cfg)
}

plot_top_forest <- function(meta_res, cfg, n = 12) {
  genes <- meta_res$meta |>
    dplyr::filter(.data$fdr < 0.05, is.finite(.data$meta_logFC)) |>
    dplyr::slice_head(n = n) |>
    dplyr::pull(.data$feature)
  if (!length(genes)) {
    return(save_status_figure("forest_gene_no_significant_genes.png", "No FDR-significant finite meta-analysis genes available for forest plots.", cfg))
  }
  purrr::map_chr(genes, plot_forest_gene, meta_res = meta_res, cfg = cfg)
}

plot_deconv_by_condition <- function(deconv_list, dataset_list, cfg) {
  safe_figure({
    df <- purrr::imap_dfr(deconv_list, function(sc, dsid) {
      prepared <- prepare_dataset_for_downstream_scores(dataset_list[[dsid]])
      meta <- prepared$meta
      sc |>
        tibble::rownames_to_column("sample") |>
        dplyr::left_join(meta |> dplyr::select("sample", "condition"), by = "sample") |>
        tidyr::pivot_longer(-c("sample", "condition"), names_to = "score", values_to = "value") |>
        dplyr::filter(.data$score %in% c("stromal_caf", "immune_pan", "ductal_epithelial", "acinar_pancreas")) |>
        dplyr::filter(is.finite(.data$value), !is.na(.data$condition)) |>
        dplyr::mutate(dataset = dsid)
    })
    if (!nrow(df)) stop("No finite deconvolution scores available for plotting.")
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$condition, .data$value, fill = .data$condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
      ggplot2::geom_jitter(width = 0.15, size = 0.6, alpha = 0.35) +
      ggplot2::facet_grid(.data$score ~ .data$dataset, scales = "free_y") +
      ggplot2::labs(x = NULL, y = "Signature score", title = "Tumor microenvironment signature scores by cohort") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = ggplot2::element_text(face = "bold"))
    save_plot(p, "deconvolution_signature_scores_by_condition.png", width = 12, height = 7, cfg = cfg)
  }, filename = "deconvolution_signature_scores_by_condition.png", cfg = cfg)
}

plot_transition_top <- function(dataset_list, transition_res, cfg, n = 20) {
  safe_figure({
    if (!nrow(transition_res) || !"GSE91035" %in% names(dataset_list)) stop("GSE91035 transition results are unavailable.")
    genes <- transition_res |>
      dplyr::filter(.data$monotonic_class == "benign_to_pdac_transition") |>
      dplyr::filter(is.finite(.data$trend_p) | is.finite(.data$pdac_vs_benign_logFC)) |>
      dplyr::slice_head(n = n) |>
      dplyr::pull(.data$gene)
    if (!length(genes)) stop("No transition genes available for heatmap.")
    ds <- dataset_list[["GSE91035"]]
    prepared <- prepare_dataset_for_downstream_scores(ds)
    ds$expr <- prepared$expr
    ds$meta <- prepared$meta
    meta <- ds$meta |> dplyr::filter(.data$condition %in% c("control", "benign", "tumor"))
    genes_present <- intersect(genes, rownames(ds$expr))
    if (length(genes_present) < 2) stop("Fewer than two transition genes are present in GSE91035 expression matrix.")
    z <- zscore_rows(ds$expr[genes_present, meta$sample, drop = FALSE])
    z <- prepare_matrix_for_heatmap(z, min_finite_per_row = 2, min_finite_per_col = 2)
    meta <- meta |> dplyr::filter(.data$sample %in% colnames(z))
    z <- z[, meta$sample, drop = FALSE]
    ann <- meta |> dplyr::select("sample", "condition") |> tibble::column_to_rownames("sample")
    path <- file.path(cfg$project$output_dir, "figures", "gse91035_transition_heatmap_top_genes.png")
    grDevices::png(path, width = 1800, height = 2200, res = 300)
    pheatmap::pheatmap(
      z,
      annotation_col = ann,
      cluster_cols = ncol(z) > 1,
      cluster_rows = nrow(z) > 1,
      show_colnames = FALSE,
      fontsize_row = 6,
      main = "Normal-benign-PDAC transition genes"
    )
    grDevices::dev.off()
    path
  }, filename = "gse91035_transition_heatmap_top_genes.png", cfg = cfg)
}

make_all_figures <- function(dataset_list, meta_res, deconv_list, transition_res, cfg) {
  fig_paths <- c(
    plot_all_qc(dataset_list, cfg),
    plot_replication_heatmap(meta_res, cfg),
    plot_heterogeneity(meta_res, cfg),
    plot_top_forest(meta_res, cfg),
    plot_deconv_by_condition(deconv_list, dataset_list, cfg),
    plot_transition_top(dataset_list, transition_res, cfg)
  )
  tibble::tibble(figure = basename(fig_paths), path = fig_paths) |>
    dplyr::filter(!is.na(.data$path)) |>
    write_tsv(file.path(cfg$project$output_dir, "tables", "generated_figures.tsv"))
  fig_paths
}
