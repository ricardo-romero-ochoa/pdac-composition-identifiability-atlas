# Tumor microenvironment scoring and confounder-aware DE ----------------------

builtin_tme_signatures <- function() {
  list(
    stromal_caf = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A3", "ACTA2", "TAGLN", "FAP", "PDGFRB", "LUM", "DCN", "THY1", "SPARC", "FN1", "MMP2", "MMP11"),
    immune_pan = c("PTPRC", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "LST1", "TYROBP", "FCER1G", "AIF1", "ITGAM", "CD3D", "CD3E", "CD68"),
    myeloid_macrophage = c("CD68", "CSF1R", "LST1", "AIF1", "TYROBP", "C1QA", "C1QB", "C1QC", "MSR1", "FCGR3A", "MRC1", "CD163"),
    t_cell = c("CD3D", "CD3E", "CD2", "TRAC", "LCK", "CD247", "CD8A", "CD8B", "IL7R", "CCR7"),
    ductal_epithelial = c("KRT19", "KRT8", "KRT18", "KRT7", "EPCAM", "MUC1", "TACSTD2", "CEACAM6", "S100P", "MMP7"),
    acinar_pancreas = c("PRSS1", "PRSS2", "CPA1", "CPA2", "CPB1", "CTRB1", "CTRB2", "CLPS", "CELA3A", "CELA3B", "REG1A", "PNLIP"),
    endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "ENG", "ESAM", "CDH5", "RAMP2"),
    proliferation = c("MKI67", "TOP2A", "BUB1", "BUB1B", "CDK1", "CCNB1", "CCNB2", "AURKA", "AURKB", "MCM2", "MCM5")
  )
}

score_signatures <- function(expr, signatures = builtin_tme_signatures()) {
  expr <- as_numeric_matrix(expr)
  if (is.null(colnames(expr)) || anyNA(colnames(expr)) || any(!nzchar(colnames(expr)))) {
    stop("Expression matrix passed to score_signatures() must have non-empty sample column names.", call. = FALSE)
  }
  if (anyDuplicated(colnames(expr))) {
    dup <- unique(colnames(expr)[duplicated(colnames(expr))])
    stop(
      "Expression matrix passed to score_signatures() has duplicated sample names: ",
      paste(utils::head(dup, 10), collapse = ", "),
      if (length(dup) > 10) " ..." else "",
      ". Average technical replicates or make sample IDs unique before scoring.",
      call. = FALSE
    )
  }

  z <- zscore_rows(expr)
  scores_long <- purrr::map_dfr(names(signatures), function(nm) {
    genes <- intersect(signatures[[nm]], rownames(z))
    if (length(genes) < 3) {
      vals <- rep(NA_real_, ncol(z))
    } else {
      vals <- colMeans(z[genes, , drop = FALSE], na.rm = TRUE)
    }
    tibble::tibble(signature = nm, sample = colnames(z), score = vals, n_genes = length(genes))
  })

  # IMPORTANT: n_genes differs by signature. If pivot_wider() is allowed to
  # infer id columns, it treats both sample and n_genes as identifiers, creating
  # multiple rows per sample and duplicate row.names. Keep n_genes as an
  # attribute/side table and widen strictly by sample.
  n_genes <- scores_long |>
    dplyr::distinct(.data$signature, .data$n_genes)
  scores <- scores_long |>
    dplyr::select("sample", "signature", "score") |>
    tidyr::pivot_wider(id_cols = "sample", names_from = "signature", values_from = "score") |>
    tibble::column_to_rownames("sample") |>
    as.data.frame()
  attr(scores, "n_genes") <- n_genes
  scores
}

run_xcell_optional <- function(expr) {
  if (!quiet_require("xCell")) return(NULL)
  out <- tryCatch({
    xCell::xCellAnalysis(expr)
  }, error = function(e) {
    warning("xCell failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(out)) return(NULL)
  t(out) |> as.data.frame()
}

prepare_dataset_for_downstream_scores <- function(ds) {
  # Use exactly the same expression/metadata alignment as limma. This matters
  # for GSE15471, where technical replicate arrays are averaged before modeling.
  average_technical_replicates(
    ds$expr,
    ds$meta,
    enabled = isTRUE(ds$config$technical_replicates)
  )
}

run_deconvolution_one <- function(ds) {
  prepared <- prepare_dataset_for_downstream_scores(ds)
  basic <- score_signatures(prepared$expr)
  xcell <- run_xcell_optional(prepared$expr)
  if (!is.null(xcell)) {
    common <- intersect(rownames(basic), rownames(xcell))
    basic <- cbind(basic[common, , drop = FALSE], xcell[common, , drop = FALSE])
  }
  rownames(basic) <- make.unique(rownames(basic))
  basic
}

run_all_deconvolution <- function(dataset_list, cfg) {
  out <- purrr::map(dataset_list, run_deconvolution_one)
  long <- purrr::imap_dfr(out, function(score, dsid) {
    score |>
      tibble::rownames_to_column("sample") |>
      tidyr::pivot_longer(-"sample", names_to = "score", values_to = "value") |>
      dplyr::mutate(dataset = dsid, .before = 1)
  })
  write_tsv(long, file.path(cfg$project$output_dir, "tables", "deconvolution_scores_long.tsv"))
  out
}

tme_adjustment_specs <- function(cfg = NULL) {
  # v3.14: pre-declared sensitivity grid. S4 is the historical, deliberately
  # over-adjusted specification and should be interpreted as a negative-control
  # adjustment because epithelial/acinar covariates are disease-state proxies in
  # tumor-versus-control contrasts.
  list(
    S0 = character(),
    S1 = c("stromal_caf", "immune_pan"),
    S2 = c("stromal_caf", "immune_pan", "endothelial"),
    S3 = "COMPOSITION_PC1_NON_EPITHELIAL",
    S4 = c("stromal_caf", "immune_pan", "ductal_epithelial", "acinar_pancreas")
  )
}

tme_spec_table <- function(cfg = NULL) {
  specs <- tme_adjustment_specs(cfg)
  tibble::tibble(
    tme_spec = names(specs),
    covariates_requested = vapply(specs, paste, collapse = ";", FUN.VALUE = character(1)),
    interpretation = c(
      "unadjusted_reference",
      "stromal_immune_only_no_epithelial_proxies",
      "stromal_immune_endothelial_no_epithelial_proxies",
      "single_non_epithelial_composition_scalar",
      "historical_over_adjusted_negative_control_with_epithelial_acinar_proxies"
    )
  )
}

make_tme_covariates_for_spec <- function(sc, spec_id, sample_order, cfg = NULL) {
  sc <- as.data.frame(sc)
  sc <- sc[sample_order, , drop = FALSE]
  specs <- tme_adjustment_specs(cfg)
  if (!spec_id %in% names(specs)) stop("Unknown TME adjustment spec: ", spec_id, call. = FALSE)
  requested <- specs[[spec_id]]
  if (!length(requested)) {
    return(NULL)
  }
  if (identical(requested, "COMPOSITION_PC1_NON_EPITHELIAL")) {
    pc_cols <- intersect(c("stromal_caf", "immune_pan", "myeloid_macrophage", "t_cell", "endothelial"), colnames(sc))
    pc_mat <- sc[, pc_cols, drop = FALSE]
    pc_mat <- pc_mat[, colSums(is.finite(as.matrix(pc_mat))) >= 5, drop = FALSE]
    if (ncol(pc_mat) < 2) {
      return(NULL)
    }
    pc_mat <- as.data.frame(scale(pc_mat))
    pc_mat[!is.finite(as.matrix(pc_mat))] <- NA_real_
    ok <- stats::complete.cases(pc_mat)
    vals <- rep(NA_real_, nrow(pc_mat))
    if (sum(ok) >= 5) {
      pr <- stats::prcomp(pc_mat[ok, , drop = FALSE], center = FALSE, scale. = FALSE)
      vals[ok] <- as.numeric(pr$x[, 1])
      # Orient the scalar so higher values tend to mean more non-epithelial
      # admixture, using stromal_caf when present as an anchor.
      if ("stromal_caf" %in% colnames(sc)) {
        anchor <- sc[["stromal_caf"]]
        cc <- suppressWarnings(stats::cor(vals, anchor, method = "spearman", use = "pairwise.complete.obs"))
        if (is.finite(cc) && cc < 0) vals <- -vals
      }
    }
    out <- data.frame(composition_pc1_non_epithelial = vals, row.names = rownames(sc), check.names = FALSE)
    attr(out, "covariates_used") <- pc_cols
    return(out)
  }
  covar_names <- intersect(requested, colnames(sc))
  if (!length(covar_names)) {
    return(NULL)
  }
  out <- sc[, covar_names, drop = FALSE]
  attr(out, "covariates_used") <- covar_names
  out
}

run_tme_adjusted_de_for_spec <- function(dataset_list, deconv_list, cfg, spec_id) {
  purrr::imap(dataset_list, function(ds, dsid) {
    sc <- deconv_list[[dsid]]
    prepared <- prepare_dataset_for_downstream_scores(ds)
    missing <- setdiff(prepared$meta$sample, rownames(sc))
    if (length(missing)) {
      stop(
        dsid, ": deconvolution scores are not aligned to prepared expression samples. Missing scores for: ",
        paste(utils::head(missing, 10), collapse = ", "),
        if (length(missing) > 10) " ..." else "",
        ". This usually indicates technical-replicate averaging was not applied consistently.",
        call. = FALSE
      )
    }
    covars <- make_tme_covariates_for_spec(sc, spec_id, prepared$meta$sample, cfg)
    used <- attr(covars, "covariates_used") %||% character()
    if (spec_id != "S0" && (is.null(covars) || !ncol(covars))) {
      warning(dsid, ": no covariates available for ", spec_id, "; returning unadjusted DE for this spec.")
      covars <- NULL
    }
    run_de_one(ds, cfg, contrast = "tumor_vs_control", adjusted_scores = covars) |>
      dplyr::mutate(
        tme_spec = spec_id,
        tme_covariates_used = paste(used, collapse = ";"),
        .after = "platform"
      )
  })
}

run_stroma_adjusted_de <- function(dataset_list, deconv_list, cfg) {
  # Backward-compatible target used by downstream code. In v3.14 the primary
  # adjusted branch is S1 (stromal/immune only). The historical S4 branch is
  # still produced by run_tme_adjustment_grid() and interpreted as an
  # over-adjusted negative control.
  primary_spec <- cfg$analysis$tme_adjustment$primary_spec %||% "S1"
  out <- run_tme_adjusted_de_for_spec(dataset_list, deconv_list, cfg, primary_spec)
  all <- dplyr::bind_rows(out)
  write_tsv(all, file.path(cfg$project$output_dir, "tables", "within_study_de_tme_adjusted.tsv"))
  write_tsv(all, file.path(cfg$project$output_dir, "tables", paste0("within_study_de_tme_adjusted_", primary_spec, ".tsv")))
  out
}

run_tme_adjustment_grid <- function(dataset_list, deconv_list, cfg) {
  specs <- tme_adjustment_specs(cfg)
  spec_info <- tme_spec_table(cfg)
  write_tsv(spec_info, file.path(cfg$project$output_dir, "tables", "tme_adjustment_specifications.tsv"))

  per_spec <- purrr::map(names(specs), function(spec_id) run_tme_adjusted_de_for_spec(dataset_list, deconv_list, cfg, spec_id))
  names(per_spec) <- names(specs)
  all <- purrr::imap_dfr(per_spec, ~ dplyr::bind_rows(.x))
  write_tsv(all, file.path(cfg$project$output_dir, "tables", "within_study_de_tme_adjustment_grid.tsv"))

  method <- cfg$analysis$meta_method %||% "REML"
  meta <- all |>
    dplyr::group_by(.data$tme_spec, .data$gene) |>
    dplyr::group_modify(~ meta_one_feature(.x, method = method)) |>
    dplyr::ungroup() |>
    dplyr::rename(feature = .data$gene) |>
    dplyr::group_by(.data$tme_spec) |>
    dplyr::mutate(fdr = stats::p.adjust(.data$p, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::left_join(spec_info, by = "tme_spec") |>
    dplyr::arrange(.data$tme_spec, .data$fdr, dplyr::desc(abs(.data$meta_logFC)))
  write_tsv(meta, file.path(cfg$project$output_dir, "tables", "gene_tme_adjustment_grid_meta.tsv"))

  base <- meta |>
    dplyr::filter(.data$tme_spec == "S0") |>
    dplyr::select(feature, unadjusted_meta_logFC = .data$meta_logFC, unadjusted_fdr = .data$fdr)
  retained_fdr <- cfg$analysis$tme_classification$retained_adjusted_fdr %||% 0.05
  unadj_fdr <- cfg$analysis$tme_classification$unadjusted_fdr %||% 0.05
  classified <- meta |>
    dplyr::left_join(base, by = "feature") |>
    dplyr::mutate(
      same_direction_as_unadjusted = is.finite(.data$meta_logFC) & is.finite(.data$unadjusted_meta_logFC) & sign(.data$meta_logFC) == sign(.data$unadjusted_meta_logFC),
      attenuation_vs_unadjusted = 1 - abs(.data$meta_logFC) / pmax(abs(.data$unadjusted_meta_logFC), 1e-8),
      attenuation_vs_unadjusted = pmax(pmin(.data$attenuation_vs_unadjusted, 1), -1),
      retained_under_spec = .data$unadjusted_fdr < unadj_fdr & .data$fdr < retained_fdr & .data$same_direction_as_unadjusted,
      attenuated_ge_075 = .data$unadjusted_fdr < unadj_fdr & is.finite(.data$attenuation_vs_unadjusted) & .data$attenuation_vs_unadjusted >= 0.75
    )
  write_tsv(classified, file.path(cfg$project$output_dir, "tables", "gene_tme_adjustment_grid_classified.tsv"))

  summary <- classified |>
    dplyr::group_by(.data$tme_spec, .data$covariates_requested, .data$interpretation) |>
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_unadjusted_meta_significant = sum(.data$unadjusted_fdr < unadj_fdr, na.rm = TRUE),
      n_retained_under_spec = sum(.data$retained_under_spec, na.rm = TRUE),
      n_attenuated_ge_075 = sum(.data$attenuated_ge_075, na.rm = TRUE),
      median_attenuation = stats::median(.data$attenuation_vs_unadjusted[.data$unadjusted_fdr < unadj_fdr], na.rm = TRUE),
      .groups = "drop"
    )
  s1_n <- summary$n_retained_under_spec[summary$tme_spec == "S1"]
  if (!length(s1_n) || !is.finite(s1_n) || s1_n == 0) s1_n <- NA_real_
  summary <- summary |>
    dplyr::mutate(retained_fold_vs_S1 = .data$n_retained_under_spec / s1_n)
  write_tsv(summary, file.path(cfg$project$output_dir, "tables", "tme_adjustment_specification_summary.tsv"))

  list(per_spec = per_spec, all = all, meta = meta, classified = classified, summary = summary, specifications = spec_info)
}

gene_deconv_cor_one <- function(ds, deconv_scores, cfg = NULL) {
  prepared <- prepare_dataset_for_downstream_scores(ds)
  expr <- prepared$expr
  meta <- prepared$meta
  missing <- setdiff(meta$sample, rownames(deconv_scores))
  if (length(missing)) {
    stop(
      ds$dataset, ": deconvolution score rows do not match prepared metadata samples for gene-score correlations. Missing: ",
      paste(utils::head(missing, 10), collapse = ", "),
      if (length(missing) > 10) " ..." else "",
      call. = FALSE
    )
  }
  scores <- deconv_scores[meta$sample, , drop = FALSE]
  scores <- scores[, colSums(is.finite(as.matrix(scores))) >= 5, drop = FALSE]
  if (!ncol(scores)) return(tibble::tibble())
  min_n <- cfg$analysis$gene_deconv_correlation$min_samples_per_group %||% 6
  groups <- c("tumor", "control", "pooled")
  purrr::map_dfr(groups, function(group_name) {
    idx <- if (group_name == "pooled") seq_len(nrow(meta)) else which(meta$condition == group_name)
    if (length(idx) < min_n) return(tibble::tibble())
    samples <- meta$sample[idx]
    purrr::map_dfr(colnames(scores), function(score_name) {
      vals <- scores[samples, score_name, drop = TRUE]
      if (sum(is.finite(vals)) < min_n) return(tibble::tibble())
      cors <- apply(expr[, samples, drop = FALSE], 1, function(g) suppressWarnings(stats::cor(g, vals, method = "spearman", use = "pairwise.complete.obs")))
      tibble::tibble(gene = names(cors), score = score_name, condition_group = group_name, n_samples = length(samples), rho = as.numeric(cors))
    })
  }) |>
    dplyr::mutate(dataset = ds$dataset, .before = 1)
}

run_gene_deconv_correlations <- function(dataset_list, deconv_list, cfg) {
  all <- purrr::imap_dfr(dataset_list, ~ gene_deconv_cor_one(.x, deconv_list[[.y]], cfg = cfg))
  write_tsv(all, file.path(cfg$project$output_dir, "tables", "gene_deconvolution_spearman_within_group.tsv"))
  # Backward-compatible filename now contains an explicit condition_group column;
  # downstream summaries use the tumor-only rows by default.
  write_tsv(all, file.path(cfg$project$output_dir, "tables", "gene_deconvolution_spearman.tsv"))
  corr_group <- cfg$analysis$gene_deconv_correlation$classification_group %||% "tumor"
  summary <- all |>
    dplyr::filter(.data$condition_group == corr_group, .data$score %in% c("stromal_caf", "immune_pan")) |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(
      mean_abs_stroma_cor = mean(abs(.data$rho[.data$score == "stromal_caf"]), na.rm = TRUE),
      mean_abs_immune_cor = mean(abs(.data$rho[.data$score == "immune_pan"]), na.rm = TRUE),
      n_stroma_cor_datasets = sum(.data$score == "stromal_caf" & is.finite(.data$rho)),
      n_immune_cor_datasets = sum(.data$score == "immune_pan" & is.finite(.data$rho)),
      correlation_group = corr_group,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      mean_abs_stroma_cor = dplyr::if_else(is.nan(.data$mean_abs_stroma_cor), NA_real_, .data$mean_abs_stroma_cor),
      mean_abs_immune_cor = dplyr::if_else(is.nan(.data$mean_abs_immune_cor), NA_real_, .data$mean_abs_immune_cor)
    )
  write_tsv(summary, file.path(cfg$project$output_dir, "tables", "gene_deconvolution_correlation_summary.tsv"))

  audit <- all |>
    dplyr::filter(.data$score %in% c("stromal_caf", "immune_pan")) |>
    dplyr::group_by(.data$condition_group, .data$score) |>
    dplyr::summarise(n_gene_score_correlations = sum(is.finite(.data$rho)), median_abs_rho = stats::median(abs(.data$rho), na.rm = TRUE), .groups = "drop")
  write_tsv(audit, file.path(cfg$project$output_dir, "tables", "gene_deconvolution_correlation_group_audit.tsv"))
  summary
}
