# Manuscript reinforcement layer ---------------------------------------------
# This file turns the exploratory atlas into a more defensible, manuscript-ready
# analysis by auditing metadata, tiering signatures, separating tumor-intrinsic
# and TME-sensitive biology, validating the GSE91035 transition axis externally,
# and re-ranking hubs with explicit penalties for weak cross-cohort support.

metadata_expected_counts <- function(cfg) {
  purrr::imap_dfr(cfg$datasets, function(ds_cfg, dsid) {
    exp <- ds_cfg$expected_counts %||% NULL
    if (is.null(exp)) return(tibble::tibble())
    tibble::tibble(
      dataset = dsid,
      condition = names(exp),
      expected_n = as.integer(unlist(exp, use.names = FALSE))
    )
  })
}

metadata_audit_report <- function(metadata_list, cfg) {
  meta <- dplyr::bind_rows(metadata_list)
  out_dir <- file.path(cfg$project$output_dir, "tables")

  observed <- meta |>
    dplyr::mutate(condition = dplyr::coalesce(.data$condition, "unlabeled")) |>
    dplyr::count(.data$dataset, .data$condition, name = "observed_n") |>
    dplyr::arrange(.data$dataset, .data$condition)

  expected <- metadata_expected_counts(cfg)
  if (nrow(expected)) {
    comparison <- dplyr::full_join(expected, observed, by = c("dataset", "condition")) |>
      dplyr::mutate(
        expected_n = dplyr::coalesce(.data$expected_n, 0L),
        observed_n = dplyr::coalesce(.data$observed_n, 0L),
        delta = .data$observed_n - .data$expected_n,
        status = dplyr::case_when(
          .data$condition == "unlabeled" & .data$observed_n > 0 ~ "review_unlabeled_samples",
          .data$delta == 0 ~ "ok",
          TRUE ~ "review_count_mismatch"
        )
      ) |>
      dplyr::arrange(.data$dataset, .data$condition)
  } else {
    comparison <- observed |>
      dplyr::mutate(expected_n = NA_integer_, delta = NA_integer_, status = "no_expected_counts_in_config") |>
      dplyr::select(.data$dataset, .data$condition, .data$expected_n, .data$observed_n, .data$delta, .data$status)
  }

  samples_to_review <- meta |>
    dplyr::mutate(
      review_reason = dplyr::case_when(
        is.na(.data$condition) | !nzchar(.data$condition) ~ "missing_condition",
        is.na(.data$include) | !.data$include ~ "excluded_by_rule_or_override",
        .data$condition %in% c("benign", "control", "tumor", "cp") ~ NA_character_,
        TRUE ~ "unexpected_condition_label"
      )
    ) |>
    dplyr::filter(!is.na(.data$review_reason)) |>
    dplyr::select(dplyr::any_of(c("dataset", "sample", "title", "condition", "include", "patient_id", "technical_group", "notes", "review_reason")))

  write_tsv(observed, file.path(out_dir, "metadata_observed_counts.tsv"))
  write_tsv(comparison, file.path(out_dir, "metadata_expected_vs_observed.tsv"))
  write_tsv(samples_to_review, file.path(out_dir, "metadata_samples_requiring_review.tsv"))

  fail_on_mismatch <- isTRUE(cfg$analysis$metadata_audit$fail_on_mismatch %||% FALSE)
  if (fail_on_mismatch && any(comparison$status != "ok")) {
    bad <- comparison |> dplyr::filter(.data$status != "ok")
    stop(
      "Metadata audit failed. Review results/tables/metadata_expected_vs_observed.tsv. First mismatch: ",
      paste(utils::capture.output(print(utils::head(bad, 5))), collapse = " | "),
      call. = FALSE
    )
  }

  list(observed = observed, comparison = comparison, samples_to_review = samples_to_review)
}

classify_signature_tiers <- function(meta_res, cfg) {
  rules <- cfg$analysis$signature_tiers %||% list()
  tier1 <- rules$tier1 %||% list()
  tier2 <- rules$tier2 %||% list()

  meta <- meta_res$meta |>
    dplyr::mutate(
      lodo_direction_concordance = dplyr::coalesce(.data$lodo_direction_concordance, .data$direction_concordance),
      lodo_max_fdr_proxy_p = dplyr::coalesce(.data$lodo_max_fdr_proxy_p, .data$p),
      signature_tier = dplyr::case_when(
        .data$fdr <= (tier1$fdr %||% 0.05) &
          abs(.data$meta_logFC) >= (tier1$abs_logfc %||% 0.40) &
          .data$k >= (tier1$min_k %||% 5) &
          .data$direction_concordance >= (tier1$min_direction_concordance %||% 1.0) &
          .data$lodo_direction_concordance >= (tier1$min_lodo_direction_concordance %||% 1.0) &
          .data$lodo_max_fdr_proxy_p <= (tier1$max_lodo_proxy_p %||% 0.10) &
          (is.na(.data$I2) | .data$I2 <= (tier1$max_i2 %||% 75)) ~ "tier1_core_cross_cohort",
        .data$fdr <= (tier2$fdr %||% 0.05) &
          abs(.data$meta_logFC) >= (tier2$abs_logfc %||% 0.40) &
          .data$k >= (tier2$min_k %||% 3) &
          .data$direction_concordance >= (tier2$min_direction_concordance %||% 1.0) ~ "tier2_recurrent",
        .data$fdr <= (tier2$fdr %||% 0.05) & .data$k < (tier2$min_k %||% 3) ~ "platform_limited_or_dataset_specific",
        .data$fdr <= 0.05 ~ "significant_but_not_core",
        TRUE ~ "not_significant"
      ),
      program = dplyr::case_when(
        .data$signature_tier == "tier1_core_cross_cohort" & .data$meta_logFC > 0 ~ "tier1_core_up",
        .data$signature_tier == "tier1_core_cross_cohort" & .data$meta_logFC < 0 ~ "tier1_core_down",
        .data$signature_tier == "tier2_recurrent" & .data$meta_logFC > 0 ~ "tier2_recurrent_up",
        .data$signature_tier == "tier2_recurrent" & .data$meta_logFC < 0 ~ "tier2_recurrent_down",
        TRUE ~ .data$signature_tier
      ),
      main_text_signature = .data$signature_tier == "tier1_core_cross_cohort"
    ) |>
    dplyr::arrange(factor(.data$signature_tier, levels = c("tier1_core_cross_cohort", "tier2_recurrent", "platform_limited_or_dataset_specific", "significant_but_not_core", "not_significant")), .data$fdr)

  write_tsv(meta, file.path(cfg$project$output_dir, "tables", "tiered_core_signature.tsv"))
  write_tsv(meta |> dplyr::filter(.data$signature_tier == "tier1_core_cross_cohort"), file.path(cfg$project$output_dir, "tables", "program_tier1_core_signature.tsv"))
  write_tsv(meta |> dplyr::filter(.data$signature_tier == "tier2_recurrent"), file.path(cfg$project$output_dir, "tables", "program_tier2_recurrent_signature.tsv"))
  write_tsv(meta |> dplyr::filter(.data$signature_tier == "platform_limited_or_dataset_specific"), file.path(cfg$project$output_dir, "tables", "program_platform_limited_signature.tsv"))
  meta
}

tme_retention_classification <- function(meta_res, adjusted_meta, gene_deconv_cor, cfg) {
  rules <- cfg$analysis$tme_classification %||% list()
  adj <- adjusted_meta$meta |>
    dplyr::select(feature, adjusted_meta_logFC = .data$meta_logFC, adjusted_fdr = .data$fdr, adjusted_k = .data$k)
  cor <- gene_deconv_cor |>
    dplyr::select(feature = .data$gene, .data$mean_abs_stroma_cor, .data$mean_abs_immune_cor)

  out <- meta_res$meta |>
    dplyr::left_join(adj, by = "feature") |>
    dplyr::left_join(cor, by = "feature") |>
    dplyr::mutate(
      adjusted_same_direction = !is.na(.data$adjusted_meta_logFC) & sign(.data$adjusted_meta_logFC) == sign(.data$meta_logFC),
      retained_after_tme_adjustment = !is.na(.data$adjusted_fdr) & .data$adjusted_fdr < (rules$retained_adjusted_fdr %||% 0.05) & .data$adjusted_same_direction,
      stromal_attenuation = 1 - abs(.data$adjusted_meta_logFC) / pmax(abs(.data$meta_logFC), 1e-8),
      stromal_attenuation = pmax(pmin(.data$stromal_attenuation, 1), -1),
      mean_abs_stroma_cor = dplyr::coalesce(.data$mean_abs_stroma_cor, 0),
      mean_abs_immune_cor = dplyr::coalesce(.data$mean_abs_immune_cor, 0),
      tme_class = dplyr::case_when(
        .data$fdr >= (rules$unadjusted_fdr %||% 0.05) ~ "not_meta_significant",
        .data$retained_after_tme_adjustment & .data$stromal_attenuation < (rules$tumor_intrinsic_max_attenuation %||% 0.75) ~ "tumor_associated_retained_under_S1",
        !.data$retained_after_tme_adjustment & (
          .data$stromal_attenuation >= (rules$microenvironment_min_attenuation %||% 0.75) |
            .data$mean_abs_stroma_cor >= (rules$microenvironment_min_abs_cor %||% 0.35) |
            .data$mean_abs_immune_cor >= (rules$microenvironment_min_abs_cor %||% 0.35)
        ) ~ "composition_covarying",
        TRUE ~ "mixed_or_ambiguous"
      ),
      tme_class_legacy = dplyr::case_when(
        .data$tme_class == "tumor_associated_retained_under_S1" ~ "tumor_intrinsic_retained",
        .data$tme_class == "composition_covarying" ~ "microenvironment_sensitive",
        TRUE ~ .data$tme_class
      ),
      classification_note = dplyr::case_when(
        .data$tme_class == "tumor_associated_retained_under_S1" ~ "Retained after v3.14 S1 stromal/immune-only adjustment; not proof of cell-intrinsic biology.",
        .data$tme_class == "composition_covarying" ~ "Tumor-associated gene whose effect is attenuated by S1 adjustment and/or covaries with stromal/immune scores within tumors.",
        TRUE ~ NA_character_
      ),
      program = dplyr::case_when(
        .data$tme_class == "tumor_associated_retained_under_S1" & .data$meta_logFC > 0 ~ "tumor_associated_retained_up",
        .data$tme_class == "tumor_associated_retained_under_S1" & .data$meta_logFC < 0 ~ "tumor_associated_retained_down",
        .data$tme_class == "composition_covarying" & .data$meta_logFC > 0 ~ "composition_covarying_up",
        .data$tme_class == "composition_covarying" & .data$meta_logFC < 0 ~ "composition_covarying_down",
        TRUE ~ .data$tme_class
      )
    ) |>
    dplyr::arrange(factor(.data$tme_class, levels = c("tumor_associated_retained_under_S1", "composition_covarying", "mixed_or_ambiguous", "not_meta_significant")), .data$fdr)

  write_tsv(out, file.path(cfg$project$output_dir, "tables", "tumor_tme_classified_programs.tsv"))
  write_tsv(out |> dplyr::filter(.data$tme_class == "tumor_associated_retained_under_S1"), file.path(cfg$project$output_dir, "tables", "program_tumor_intrinsic_retained_signature.tsv"))
  write_tsv(out |> dplyr::filter(.data$tme_class == "composition_covarying"), file.path(cfg$project$output_dir, "tables", "program_composition_covarying_signature.tsv"))
  # Backward-compatible alias for earlier external-validation scripts; do not
  # use the legacy name in the manuscript text.
  write_tsv(out |> dplyr::filter(.data$tme_class == "composition_covarying"), file.path(cfg$project$output_dir, "tables", "program_microenvironment_sensitive_signature.tsv"))
  write_tsv(out |> dplyr::filter(.data$tme_class == "mixed_or_ambiguous"), file.path(cfg$project$output_dir, "tables", "program_mixed_or_ambiguous_signature.tsv"))
  out
}

score_gene_set_vector <- function(expr, genes) {
  expr <- as_numeric_matrix(expr)
  genes <- intersect(unique(toupper(genes)), rownames(expr))
  if (length(genes) < 3) return(rep(NA_real_, ncol(expr)) |> stats::setNames(colnames(expr)))
  z <- zscore_rows(expr[genes, , drop = FALSE])
  out <- colMeans(z, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

module_score_test <- function(score, meta, dataset, module, expected_direction) {
  df <- tibble::tibble(sample = names(score), score = as.numeric(score)) |>
    dplyr::left_join(meta |> dplyr::select("sample", "condition"), by = "sample") |>
    dplyr::filter(.data$condition %in% c("control", "tumor"), is.finite(.data$score))
  if (dplyr::n_distinct(df$condition) < 2 || nrow(df) < 6) {
    return(tibble::tibble(dataset = dataset, module = module, n_control = sum(df$condition == "control"), n_tumor = sum(df$condition == "tumor"), mean_control = NA_real_, mean_tumor = NA_real_, logFC = NA_real_, p = NA_real_, expected_direction = expected_direction, directionally_concordant = NA))
  }
  mean_control <- mean(df$score[df$condition == "control"], na.rm = TRUE)
  mean_tumor <- mean(df$score[df$condition == "tumor"], na.rm = TRUE)
  p <- tryCatch(stats::wilcox.test(score ~ condition, data = df)$p.value, error = function(e) NA_real_)
  logFC <- mean_tumor - mean_control
  tibble::tibble(
    dataset = dataset,
    module = module,
    n_control = sum(df$condition == "control"),
    n_tumor = sum(df$condition == "tumor"),
    mean_control = mean_control,
    mean_tumor = mean_tumor,
    logFC = logFC,
    p = p,
    expected_direction = expected_direction,
    directionally_concordant = sign(logFC) == sign(expected_direction)
  )
}

validate_transition_program_external <- function(transition_signature, dataset_list, cfg) {
  if (!nrow(transition_signature)) {
    out <- tibble::tibble()
    write_tsv(out, file.path(cfg$project$output_dir, "tables", "transition_module_external_validation.tsv"))
    write_tsv(out, file.path(cfg$project$output_dir, "tables", "transition_module_external_validation_summary.tsv"))
    return(list(per_dataset = out, summary = out))
  }

  up_genes <- transition_signature |> dplyr::filter(.data$program == "transition_up") |> dplyr::pull(.data$feature) |> unique()
  down_genes <- transition_signature |> dplyr::filter(.data$program == "transition_down") |> dplyr::pull(.data$feature) |> unique()

  per_dataset <- purrr::imap_dfr(dataset_list, function(ds, dsid) {
    prepared <- prepare_dataset_for_downstream_scores(ds)
    expr <- prepared$expr
    meta <- prepared$meta
    up <- score_gene_set_vector(expr, up_genes)
    down <- score_gene_set_vector(expr, down_genes)
    composite <- up - down
    dplyr::bind_rows(
      module_score_test(up, meta, dsid, "transition_up_score", expected_direction = 1),
      module_score_test(down, meta, dsid, "transition_down_score", expected_direction = -1),
      module_score_test(composite, meta, dsid, "transition_composite_up_minus_down", expected_direction = 1)
    ) |>
      dplyr::mutate(external_to_discovery = .data$dataset != "GSE91035")
  }) |>
    dplyr::mutate(fdr = stats::p.adjust(.data$p, method = "BH"))

  validation_rules <- cfg$analysis$transition_validation %||% list()
  summary <- per_dataset |>
    dplyr::filter(.data$external_to_discovery) |>
    dplyr::group_by(.data$module) |>
    dplyr::summarise(
      n_external_tested = sum(is.finite(.data$p)),
      n_directionally_concordant = sum(.data$directionally_concordant %in% TRUE, na.rm = TRUE),
      n_nominal_significant = sum(.data$p < 0.05 & .data$directionally_concordant %in% TRUE, na.rm = TRUE),
      min_p = suppressWarnings(min(.data$p, na.rm = TRUE)),
      median_logFC = stats::median(.data$logFC, na.rm = TRUE),
      externally_supported = .data$n_directionally_concordant >= (validation_rules$min_external_concordant %||% 4),
      .groups = "drop"
    ) |>
    dplyr::mutate(min_p = dplyr::if_else(is.infinite(.data$min_p), NA_real_, .data$min_p))

  write_tsv(per_dataset, file.path(cfg$project$output_dir, "tables", "transition_module_external_validation.tsv"))
  write_tsv(summary, file.path(cfg$project$output_dir, "tables", "transition_module_external_validation_summary.tsv"))
  list(per_dataset = per_dataset, summary = summary)
}

load_manual_druggability_seed <- function(path = "resources/druggability_manual_seed.csv") {
  if (!file.exists(path)) {
    return(tibble::tibble(gene = character(), target_class = character(), subcellular_location = character(), tractability_manual = character(), known_context = character(), annotation_source = character()))
  }
  readr::read_csv(path, show_col_types = FALSE) |>
    janitor::clean_names() |>
    dplyr::mutate(gene = toupper(.data$gene))
}

build_druggability_annotation <- function(hub_table, cfg) {
  genes <- unique(toupper(utils::head(hub_table$gene, cfg$analysis$hub_prioritization$annotation_top_n %||% 500)))
  ot <- annotate_opentargets(genes, cfg) |>
    dplyr::mutate(gene = toupper(.data$symbol))
  manual <- load_manual_druggability_seed()
  ann <- tibble::tibble(gene = genes) |>
    dplyr::left_join(manual, by = "gene") |>
    dplyr::left_join(ot |> dplyr::select(-dplyr::any_of("symbol")), by = "gene")
  # Keep a stable schema when optional online annotations are disabled.
  for (cc in c("known_drug_count", "tractability", "top_known_drugs", "approved_name", "biotype", "open_targets_status", "tractability_manual")) {
    if (!cc %in% names(ann)) ann[[cc]] <- NA
  }
  ann <- ann |>
    dplyr::mutate(
      known_drug_count = suppressWarnings(as.numeric(.data$known_drug_count)),
      druggability_evidence = dplyr::case_when(
        !is.na(.data$known_drug_count) & .data$known_drug_count > 0 ~ "known_drug_or_clinical_ligand_in_open_targets",
        !is.na(.data$tractability) & nzchar(.data$tractability) ~ "open_targets_tractability_signal",
        !is.na(.data$tractability_manual) & nzchar(.data$tractability_manual) ~ .data$tractability_manual,
        TRUE ~ "no_druggability_evidence_in_current_annotations"
      )
    )
  write_tsv(ann, file.path(cfg$project$output_dir, "tables", "druggability_annotation.tsv"))
  ann
}

make_strict_integrative_hub_table <- function(hub_table, tiered_signature, tumor_tme_classes, transition_validation, cfg) {
  rules <- cfg$analysis$hub_prioritization %||% list()
  tier_info <- tiered_signature |>
    dplyr::select(feature, signature_tier, main_text_signature)
  tme_info <- tumor_tme_classes |>
    dplyr::transmute(
      feature = .data$feature,
      tme_class = .data$tme_class,
      tme_retained_after_adjustment = .data$retained_after_tme_adjustment,
      stromal_attenuation = .data$stromal_attenuation,
      mean_abs_stroma_cor = .data$mean_abs_stroma_cor,
      mean_abs_immune_cor = .data$mean_abs_immune_cor
    )

  tv_summary <- transition_validation$summary %||% tibble::tibble()
  transition_supported <- FALSE
  if (nrow(tv_summary)) {
    transition_supported <- any(tv_summary$module == "transition_composite_up_minus_down" & tv_summary$externally_supported %in% TRUE)
  }

  ann <- build_druggability_annotation(hub_table, cfg)
  annotation_cols <- c(
    "ensembl_gene_id", "approved_name", "biotype", "tractability",
    "known_drug_count", "top_known_drugs", "open_targets_status",
    "target_class", "subcellular_location", "tractability_manual",
    "known_context", "annotation_source", "druggability_evidence"
  )
  hub_base <- hub_table |> dplyr::select(-dplyr::any_of(annotation_cols))
  abs_z_rng <- range(abs(hub_base$z), na.rm = TRUE)
  if (!all(is.finite(abs_z_rng)) || diff(abs_z_rng) == 0) abs_z_rng <- c(0, 1)

  out <- hub_base |>
    dplyr::left_join(tier_info, by = c("feature" = "feature")) |>
    dplyr::left_join(tme_info, by = c("feature" = "feature")) |>
    dplyr::left_join(ann, by = "gene") |>
    dplyr::mutate(
      signature_tier = dplyr::coalesce(.data$signature_tier, "not_significant"),
      tme_class = dplyr::coalesce(.data$tme_class, "not_classified"),
      k_penalty = dplyr::case_when(
        .data$k >= (rules$main_text_min_k %||% 5) ~ 0,
        .data$k >= 3 ~ (rules$penalty_k_3_4 %||% 0.75),
        TRUE ~ (rules$penalty_k_lt3 %||% 2.00)
      ),
      heterogeneity_penalty_strict = dplyr::coalesce(.data$I2, 100) / 100,
      no_tme_retention_penalty = dplyr::if_else(.data$tme_retained_after_adjustment %in% TRUE | .data$tme_class == "composition_covarying", 0, rules$penalty_no_tme_retention %||% 0.50),
      annotation_bonus = dplyr::case_when(
        !is.na(.data$known_drug_count) & .data$known_drug_count > 0 ~ 0.60,
        !is.na(.data$tractability) & nzchar(.data$tractability) ~ 0.40,
        !is.na(.data$tractability_manual) & nzchar(.data$tractability_manual) & .data$tractability_manual != "no_known_druggability" ~ 0.25,
        TRUE ~ 0
      ),
      strict_hub_score =
        scales::rescale(abs(.data$z), to = c(0, 1), from = abs_z_rng) +
        0.80 * dplyr::coalesce(.data$direction_concordance, 0) +
        0.90 * as.numeric(.data$signature_tier == "tier1_core_cross_cohort") +
        0.60 * as.numeric(.data$tme_class == "tumor_associated_retained_under_S1") +
        0.40 * as.numeric(.data$tme_class == "composition_covarying") +
        0.30 * as.numeric(.data$transition_flag %in% TRUE & transition_supported) +
        0.25 * dplyr::coalesce(.data$module_score, 0) +
        .data$annotation_bonus -
        .data$k_penalty -
        0.65 * .data$heterogeneity_penalty_strict -
        .data$no_tme_retention_penalty,
      hub_interpretation = dplyr::case_when(
        .data$signature_tier == "tier1_core_cross_cohort" & .data$tme_class == "tumor_associated_retained_under_S1" ~ "main_text_tumor_associated_retained_candidate",
        .data$signature_tier == "tier1_core_cross_cohort" & .data$tme_class == "composition_covarying" ~ "main_text_composition_covarying_candidate",
        .data$signature_tier == "tier1_core_cross_cohort" ~ "main_text_bulk_tumor_associated_candidate",
        .data$signature_tier == "tier2_recurrent" ~ "supplementary_recurrent_candidate",
        .data$signature_tier == "platform_limited_or_dataset_specific" ~ "hypothesis_generating_platform_limited",
        TRUE ~ "not_prioritized"
      ),
      main_text_hub_candidate = .data$hub_interpretation %in% c(
        "main_text_tumor_associated_retained_candidate",
        "main_text_composition_covarying_candidate",
        "main_text_bulk_tumor_associated_candidate"
      ) & .data$k >= (rules$main_text_min_k %||% 5) & .data$fdr <= (rules$main_text_fdr %||% 0.05)
    ) |>
    dplyr::arrange(dplyr::desc(.data$main_text_hub_candidate), dplyr::desc(.data$strict_hub_score))

  write_tsv(out, file.path(cfg$project$output_dir, "tables", "integrative_hub_gene_prioritization_strict.tsv"))
  write_tsv(out |> dplyr::filter(.data$main_text_hub_candidate), file.path(cfg$project$output_dir, "tables", "main_text_hub_candidates.tsv"))
  out
}

write_reinforcement_tables <- function(metadata_audit, tiered_signature, tumor_tme_classes, transition_validation, strict_hub_table, cfg) {
  summary <- tibble::tibble(
    reinforcement = c(
      "metadata_audit_mismatches",
      "tier1_core_genes",
      "tier2_recurrent_genes",
      "platform_limited_genes",
      "tumor_associated_retained_under_S1_genes",
      "composition_covarying_genes",
      "transition_modules_externally_supported",
      "main_text_hub_candidates"
    ),
    value = c(
      sum(metadata_audit$comparison$status != "ok", na.rm = TRUE),
      sum(tiered_signature$signature_tier == "tier1_core_cross_cohort", na.rm = TRUE),
      sum(tiered_signature$signature_tier == "tier2_recurrent", na.rm = TRUE),
      sum(tiered_signature$signature_tier == "platform_limited_or_dataset_specific", na.rm = TRUE),
      sum(tumor_tme_classes$tme_class == "tumor_associated_retained_under_S1", na.rm = TRUE),
      sum(tumor_tme_classes$tme_class == "composition_covarying", na.rm = TRUE),
      if (nrow(transition_validation$summary)) sum(transition_validation$summary$externally_supported %in% TRUE, na.rm = TRUE) else 0,
      sum(strict_hub_table$main_text_hub_candidate %in% TRUE, na.rm = TRUE)
    )
  )
  write_tsv(summary, file.path(cfg$project$output_dir, "tables", "reinforcement_summary.tsv"))
  list(summary = summary)
}

plot_signature_tiers <- function(tiered_signature, cfg) {
  safe_figure({
    df <- tiered_signature |>
      dplyr::filter(.data$signature_tier != "not_significant") |>
      dplyr::count(.data$signature_tier, name = "n")
    if (!nrow(df)) stop("No signature tier counts available.")
    p <- ggplot2::ggplot(df, ggplot2::aes(stats::reorder(.data$signature_tier, .data$n), .data$n)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Genes", title = "Signature tiering for manuscript claims") +
      ggplot2::theme_bw(base_size = 11)
    save_plot(p, "reinforcement_signature_tier_counts.png", width = 6, height = 4, cfg = cfg)
  }, filename = "reinforcement_signature_tier_counts.png", cfg = cfg)
}

plot_tme_classes <- function(tumor_tme_classes, cfg) {
  safe_figure({
    df <- tumor_tme_classes |>
      dplyr::filter(.data$fdr < 0.05) |>
      dplyr::count(.data$tme_class, name = "n")
    if (!nrow(df)) stop("No TME class counts available.")
    p <- ggplot2::ggplot(df, ggplot2::aes(stats::reorder(.data$tme_class, .data$n), .data$n)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Meta-significant genes", title = "Tumor-intrinsic versus TME-sensitive classification") +
      ggplot2::theme_bw(base_size = 11)
    save_plot(p, "reinforcement_tme_class_counts.png", width = 6, height = 4, cfg = cfg)
  }, filename = "reinforcement_tme_class_counts.png", cfg = cfg)
}

plot_transition_validation <- function(transition_validation, cfg) {
  safe_figure({
    df <- transition_validation$per_dataset
    if (!nrow(df)) stop("No transition validation rows available.")
    df <- df |> dplyr::filter(.data$external_to_discovery, is.finite(.data$logFC))
    if (!nrow(df)) stop("No external transition validation rows available.")
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$dataset, .data$logFC)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_point(ggplot2::aes(shape = .data$directionally_concordant, size = -log10(pmax(.data$p, 1e-300)))) +
      ggplot2::facet_wrap(~ module, scales = "free_y") +
      ggplot2::labs(x = NULL, y = "Tumor-control score difference", title = "External validation of GSE91035 transition modules") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot(p, "reinforcement_transition_external_validation.png", width = 9, height = 4.8, cfg = cfg)
  }, filename = "reinforcement_transition_external_validation.png", cfg = cfg)
}

plot_strict_hubs <- function(strict_hub_table, cfg) {
  safe_figure({
    df <- strict_hub_table |>
      dplyr::filter(is.finite(.data$strict_hub_score)) |>
      dplyr::slice_head(n = 30) |>
      dplyr::mutate(label = paste0(.data$gene, "\n", .data$hub_interpretation))
    if (!nrow(df)) stop("No strict hub table rows available.")
    p <- ggplot2::ggplot(df, ggplot2::aes(stats::reorder(.data$gene, .data$strict_hub_score), .data$strict_hub_score)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Strict hub score", title = "Strict integrative hub prioritization") +
      ggplot2::theme_bw(base_size = 10)
    save_plot(p, "reinforcement_strict_hub_prioritization.png", width = 7, height = 7, cfg = cfg)
  }, filename = "reinforcement_strict_hub_prioritization.png", cfg = cfg)
}

make_reinforcement_figures <- function(tiered_signature, tumor_tme_classes, transition_validation, strict_hub_table, cfg) {
  fig_paths <- c(
    plot_signature_tiers(tiered_signature, cfg),
    plot_tme_classes(tumor_tme_classes, cfg),
    plot_transition_validation(transition_validation, cfg),
    plot_strict_hubs(strict_hub_table, cfg)
  )
  out <- tibble::tibble(figure = basename(fig_paths), path = fig_paths) |>
    dplyr::filter(!is.na(.data$path))
  old_path <- file.path(cfg$project$output_dir, "tables", "generated_figures.tsv")
  old <- if (file.exists(old_path)) readr::read_tsv(old_path, show_col_types = FALSE) else tibble::tibble()
  dplyr::bind_rows(old, out) |>
    dplyr::distinct(.data$figure, .keep_all = TRUE) |>
    write_tsv(old_path)
  fig_paths
}
