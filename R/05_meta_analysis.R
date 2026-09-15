# Random-effects meta-analysis and sensitivity --------------------------------

meta_one_feature <- function(df, method = "REML") {
  df <- df |> dplyr::filter(is.finite(.data$logFC), is.finite(.data$se), .data$se > 0)
  k <- nrow(df)
  if (k < 2) {
    return(tibble::tibble(
      k = k, meta_logFC = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_, tau2 = NA_real_, I2 = NA_real_, QEp = NA_real_,
      direction_concordance = NA_real_, datasets = paste(df$dataset, collapse = ";")
    ))
  }
  fit <- tryCatch(
    metafor::rma.uni(yi = df$logFC, sei = df$se, method = method),
    error = function(e) NULL
  )
  direction_concordance <- max(mean(df$logFC > 0), mean(df$logFC < 0), na.rm = TRUE)
  if (is.null(fit)) {
    return(tibble::tibble(
      k = k, meta_logFC = mean(df$logFC, na.rm = TRUE), se = NA_real_, z = NA_real_, p = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_, tau2 = NA_real_, I2 = NA_real_, QEp = NA_real_,
      direction_concordance = direction_concordance, datasets = paste(df$dataset, collapse = ";")
    ))
  }
  tibble::tibble(
    k = k,
    meta_logFC = as.numeric(fit$b[1]),
    se = fit$se,
    z = fit$zval,
    p = fit$pval,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    QEp = fit$QEp,
    direction_concordance = direction_concordance,
    datasets = paste(df$dataset, collapse = ";")
  )
}

leave_one_dataset_out <- function(df, method = "REML") {
  datasets <- unique(df$dataset)
  if (length(datasets) < 3) return(tibble::tibble())
  purrr::map_dfr(datasets, function(d) {
    res <- meta_one_feature(df |> dplyr::filter(.data$dataset != d), method = method)
    res |> dplyr::mutate(left_out = d)
  })
}

meta_analyze_de <- function(de_list, cfg, feature_col = "gene", out_prefix = "gene") {
  all <- dplyr::bind_rows(de_list)
  method <- cfg$analysis$meta_method %||% "REML"
  feature_sym <- rlang::sym(feature_col)
  meta <- all |>
    dplyr::group_by(!!feature_sym) |>
    dplyr::group_modify(~ meta_one_feature(.x, method = method)) |>
    dplyr::ungroup() |>
    dplyr::rename(feature = !!feature_sym) |>
    dplyr::mutate(fdr = stats::p.adjust(.data$p, method = "BH")) |>
    dplyr::arrange(.data$fdr, dplyr::desc(abs(.data$meta_logFC)))

  lodo <- all |>
    dplyr::semi_join(meta |> dplyr::filter(.data$fdr < 0.10) |> dplyr::select(feature), by = setNames("feature", feature_col)) |>
    dplyr::group_by(!!feature_sym) |>
    dplyr::group_modify(~ leave_one_dataset_out(.x, method = method)) |>
    dplyr::ungroup() |>
    dplyr::rename(feature = !!feature_sym)

  if (nrow(lodo)) {
    lodo_summary <- lodo |>
      dplyr::group_by(.data$feature) |>
      dplyr::summarise(
        lodo_min_abs_logFC = min(abs(.data$meta_logFC), na.rm = TRUE),
        lodo_direction_concordance = max(mean(.data$meta_logFC > 0, na.rm = TRUE), mean(.data$meta_logFC < 0, na.rm = TRUE)),
        lodo_max_fdr_proxy_p = max(.data$p, na.rm = TRUE),
        .groups = "drop"
      )
    meta <- meta |> dplyr::left_join(lodo_summary, by = "feature")
  } else {
    meta$lodo_min_abs_logFC <- NA_real_
    meta$lodo_direction_concordance <- NA_real_
    meta$lodo_max_fdr_proxy_p <- NA_real_
  }

  write_tsv(meta, file.path(cfg$project$output_dir, "tables", paste0(out_prefix, "_random_effects_meta.tsv")))
  write_tsv(lodo, file.path(cfg$project$output_dir, "tables", paste0(out_prefix, "_leave_one_dataset_out.tsv")))
  list(meta = meta, lodo = lodo, per_study = all)
}

define_core_signature <- function(meta_res, cfg) {
  thr <- cfg$analysis$core_signature
  meta_res$meta |>
    dplyr::filter(
      .data$fdr <= (thr$fdr %||% 0.05),
      abs(.data$meta_logFC) >= (thr$abs_logfc %||% 0.58),
      .data$I2 <= (thr$max_i2 %||% 60),
      .data$direction_concordance >= (thr$min_direction_concordance %||% 0.75),
      is.na(.data$lodo_direction_concordance) | .data$lodo_direction_concordance >= (thr$min_lodo_direction_concordance %||% 0.75)
    ) |>
    dplyr::mutate(program = dplyr::if_else(.data$meta_logFC > 0, "core_tumor_up", "core_tumor_down")) |>
    dplyr::arrange(.data$fdr, dplyr::desc(abs(.data$meta_logFC)))
}

define_microenvironment_program <- function(unadjusted_meta, adjusted_meta, deconv_gene_cor, cfg) {
  thr <- cfg$analysis$microenvironment_program
  adj <- adjusted_meta$meta |> dplyr::select(feature, adjusted_meta_logFC = meta_logFC, adjusted_fdr = fdr)
  cor <- deconv_gene_cor |> dplyr::select(feature = gene, mean_abs_stroma_cor, mean_abs_immune_cor)
  out <- unadjusted_meta$meta |>
    dplyr::left_join(adj, by = "feature") |>
    dplyr::left_join(cor, by = "feature") |>
    dplyr::mutate(
      stromal_attenuation = 1 - abs(.data$adjusted_meta_logFC) / pmax(abs(.data$meta_logFC), 1e-8),
      stromal_attenuation = pmax(pmin(.data$stromal_attenuation, 1), -1),
      mean_abs_stroma_cor = .data$mean_abs_stroma_cor %||% NA_real_,
      mean_abs_immune_cor = .data$mean_abs_immune_cor %||% NA_real_,
      program = dplyr::case_when(
        .data$meta_logFC > 0 ~ "composition_covarying_up",
        .data$meta_logFC < 0 ~ "composition_covarying_down",
        TRUE ~ "composition_covarying"
      )
    ) |>
    dplyr::filter(
      .data$fdr < 0.05,
      .data$stromal_attenuation >= (thr$stromal_attenuation_min %||% 0.35) |
        .data$mean_abs_stroma_cor >= (thr$abs_stroma_cor_min %||% 0.35) |
        .data$mean_abs_immune_cor >= (thr$abs_immune_cor_min %||% 0.35)
    ) |>
    dplyr::arrange(dplyr::desc(.data$stromal_attenuation), .data$fdr)
  out
}

write_program_tables <- function(core, microenv, transition, cfg) {
  write_tsv(core, file.path(cfg$project$output_dir, "tables", "program_core_tumor_signature.tsv"))
  write_tsv(microenv, file.path(cfg$project$output_dir, "tables", "program_composition_covarying_unfiltered_signature.tsv"))
  # Backward-compatible legacy output name. Use composition-covarying terminology in manuscripts.
  write_tsv(microenv, file.path(cfg$project$output_dir, "tables", "program_microenvironment_signature.tsv"))
  write_tsv(transition, file.path(cfg$project$output_dir, "tables", "program_benign_to_pdac_transition.tsv"))
  # v3.11 compatibility file for the external validation layer.
  write_tsv(transition, file.path(cfg$project$output_dir, "tables", "transition_program.tsv"))
  list(core = core, composition_covarying = microenv, transition = transition)
}
