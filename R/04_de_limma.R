# Within-study differential expression ---------------------------------------

align_expr_meta <- function(expr, meta, context = "") {
  # limma requires ncol(expr) == nrow(design). Most opaque qr.qty()/subscript
  # errors ultimately come from expression/metadata misalignment, duplicated
  # sample rows, or samples that were dropped silently by curation. Fail early
  # with actionable diagnostics instead.
  if (!"sample" %in% names(meta)) {
    stop("Metadata is missing a `sample` column", context, call. = FALSE)
  }
  meta <- as.data.frame(meta)
  meta <- meta[!is.na(meta$sample) & nzchar(as.character(meta$sample)), , drop = FALSE]
  meta$sample <- as.character(meta$sample)

  dup_meta <- unique(meta$sample[duplicated(meta$sample)])
  if (length(dup_meta)) {
    stop(
      "Duplicated sample IDs in metadata", context, ": ",
      paste(utils::head(dup_meta, 10), collapse = ", "),
      if (length(dup_meta) > 10) " ..." else "",
      ". Resolve before limma fitting.",
      call. = FALSE
    )
  }

  missing_expr <- setdiff(meta$sample, colnames(expr))
  if (length(missing_expr)) {
    stop(
      "Metadata samples are absent from expression matrix", context, ": ",
      paste(utils::head(missing_expr, 10), collapse = ", "),
      if (length(missing_expr) > 10) " ..." else "",
      ". Expression columns available: ", length(colnames(expr)),
      "; metadata rows: ", nrow(meta), ".",
      call. = FALSE
    )
  }

  meta <- meta[match(intersect(colnames(expr), meta$sample), meta$sample), , drop = FALSE]
  meta <- meta[match(meta$sample, colnames(expr), nomatch = 0L) > 0L, , drop = FALSE]
  expr <- expr[, meta$sample, drop = FALSE]
  if (ncol(expr) != nrow(meta)) {
    stop(
      "Expression/metadata dimension mismatch", context, ": ncol(expr)=", ncol(expr),
      ", nrow(meta)=", nrow(meta), ".",
      call. = FALSE
    )
  }
  storage.mode(expr) <- "numeric"
  list(expr = expr, meta = tibble::as_tibble(meta))
}

average_technical_replicates <- function(expr, meta, enabled = FALSE) {
  aligned <- align_expr_meta(expr, meta, context = " before technical-replicate handling")
  expr <- aligned$expr
  meta <- aligned$meta

  # Only GSE15471 is configured to contain technical replicate arrays. Earlier
  # versions averaged duplicate `technical_group` values for every dataset,
  # which can accidentally collapse biological samples in mixed paired/unpaired
  # cohorts such as GSE16515 and then trigger limma qr.qty() dimension errors.
  if (!isTRUE(enabled)) return(list(expr = expr, meta = meta))
  if (!"technical_group" %in% names(meta)) return(list(expr = expr, meta = meta))

  dup_groups <- as.character(meta$technical_group)
  dup_groups[is.na(dup_groups) | !nzchar(dup_groups)] <- meta$sample[is.na(dup_groups) | !nzchar(dup_groups)]
  if (!anyDuplicated(dup_groups)) return(list(expr = expr, meta = meta))

  expr2 <- limma::avearrays(expr, ID = dup_groups)
  meta$technical_group <- dup_groups
  meta2 <- meta |>
    dplyr::group_by(.data$technical_group) |>
    dplyr::summarise(
      dplyr::across(dplyr::any_of(c("dataset", "condition", "patient_id")), mode_or_na),
      sample = paste(.data$sample, collapse = ";"),
      title = paste(stats::na.omit(.data$title), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::mutate(sample = .data$technical_group) |>
    dplyr::arrange(match(.data$sample, colnames(expr2)))
  expr2 <- expr2[, meta2$sample, drop = FALSE]
  align_expr_meta(expr2, meta2, context = " after technical-replicate averaging")
}

condition_factor_for_contrast <- function(meta, contrast = "tumor_vs_control") {
  if (contrast %in% c("tumor_vs_control", "pdac_vs_normal")) {
    keep <- meta$condition %in% c("tumor", "control")
    y <- factor(ifelse(meta$condition == "tumor", "tumor", "control"), levels = c("control", "tumor"))
    return(list(keep = keep, y = y, reference_level = "control", case_level = "tumor"))
  }
  if (contrast == "pdac_vs_benign") {
    keep <- meta$condition %in% c("tumor", "benign")
    y <- factor(ifelse(meta$condition == "tumor", "tumor", "benign"), levels = c("benign", "tumor"))
    return(list(keep = keep, y = y, reference_level = "benign", case_level = "tumor"))
  }
  if (contrast == "benign_vs_normal") {
    keep <- meta$condition %in% c("benign", "control")
    y <- factor(ifelse(meta$condition == "benign", "benign", "control"), levels = c("control", "benign"))
    return(list(keep = keep, y = y, reference_level = "control", case_level = "benign"))
  }
  stop("Unknown contrast: ", contrast)
}

condition_column_name <- function(level) {
  make.names(paste0("condition", level))
}

validate_design_expr <- function(design, expr, meta, context = "") {
  if (nrow(design) != ncol(expr)) {
    stop(
      "Design/expression mismatch", context, ": nrow(design)=", nrow(design),
      ", ncol(expr)=", ncol(expr), ", nrow(meta)=", nrow(meta), ". ",
      "This usually means model.matrix dropped rows because of missing covariates, ",
      "or expression columns are not aligned to metadata samples.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

stop_missing_design_column <- function(required, design, contrast, context = "") {
  missing <- setdiff(required, colnames(design))
  if (length(missing)) {
    stop(
      "Could not build contrast ", contrast, context, ". Missing design column(s): ",
      paste(missing, collapse = ", "), ". Available design columns: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }
}

make_paired_contrast_matrix <- function(design, condition, contrast) {
  case_col <- condition_column_name(levels(condition)[2])
  ref_col <- condition_column_name(levels(condition)[1])
  stop_missing_design_column(c(case_col, ref_col), design, contrast, context = " for paired model")
  contrast_vec <- rep(0, ncol(design))
  names(contrast_vec) <- colnames(design)
  contrast_vec[case_col] <- 1
  contrast_vec[ref_col] <- -1
  mat <- matrix(contrast_vec, ncol = 1)
  rownames(mat) <- colnames(design)
  colnames(mat) <- contrast
  mat
}

coefficient_for_unpaired_model <- function(design, condition, contrast) {
  case_col <- condition_column_name(levels(condition)[2])
  stop_missing_design_column(case_col, design, contrast, context = " for unpaired or duplicate-correlation model")
  case_col
}

safe_top_table <- function(fit, coef_name, contrast, design = NULL) {
  available <- colnames(fit$coefficients)
  if (is.null(available)) available <- character()
  if (is.character(coef_name) && !coef_name %in% available) {
    if (ncol(fit$coefficients) == 1L) {
      coef_name <- 1L
    } else {
      msg <- paste0(
        "Coefficient '", coef_name, "' was not found after fitting contrast '", contrast, "'. ",
        "Available fitted coefficients: ", paste(available, collapse = ", "), "."
      )
      if (!is.null(design)) {
        msg <- paste0(msg, " Design columns: ", paste(colnames(design), collapse = ", "), ".")
      }
      stop(msg, call. = FALSE)
    }
  }
  limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
}

ordinary_se_from_limma <- function(fit_unmoderated, coef_name) {
  # v3.14: random-effects meta-analysis must not use eBayes-moderated
  # statistics to reconstruct standard errors. limma stores the ordinary
  # (pre-eBayes) residual sigma and coefficient unscaled SEs in the fitted
  # object. For contrast fits, call this on the object returned by
  # contrasts.fit() before eBayes().
  if (is.null(fit_unmoderated$stdev.unscaled) || is.null(fit_unmoderated$sigma)) {
    return(rep(NA_real_, nrow(fit_unmoderated$coefficients)) |>
      stats::setNames(rownames(fit_unmoderated$coefficients)))
  }
  coef_index <- coef_name
  if (is.character(coef_name)) {
    cn <- colnames(fit_unmoderated$coefficients)
    if (!is.null(cn) && coef_name %in% cn) coef_index <- match(coef_name, cn) else coef_index <- 1L
  }
  su <- fit_unmoderated$stdev.unscaled[, coef_index]
  se <- abs(as.numeric(su) * as.numeric(fit_unmoderated$sigma))
  se[!is.finite(se) | se <= 0] <- NA_real_
  stats::setNames(se, rownames(fit_unmoderated$coefficients))
}

run_limma_contrast <- function(expr, meta, contrast = "tumor_vs_control", paired_design = FALSE, adjust_covariates = NULL) {
  aligned0 <- align_expr_meta(expr, meta, context = paste0(" for ", contrast, " input"))
  expr <- aligned0$expr
  meta <- aligned0$meta

  cf <- condition_factor_for_contrast(meta, contrast)
  meta <- meta[cf$keep, , drop = FALSE]
  aligned <- align_expr_meta(expr, meta, context = paste0(" for ", contrast, " after condition filtering"))
  expr <- aligned$expr
  meta <- aligned$meta
  # Build the binary model condition explicitly and force the modelling
  # column to use the requested reference/case order. Earlier versions kept
  # GEO-derived `meta$condition` as a character column and created a separate
  # local variable named `condition`. In formulas such as `~ condition`,
  # model.matrix() preferentially used the metadata column, not the local
  # factor, so contrasts involving benign/control could silently get the wrong
  # baseline (e.g. conditioncontrol instead of conditionbenign).
  meta$condition_raw <- meta$condition
  condition <- droplevels(factor(
    ifelse(meta$condition_raw == cf$case_level, cf$case_level, cf$reference_level),
    levels = c(cf$reference_level, cf$case_level)
  ))
  meta$condition <- condition
  meta$condition_binary <- condition

  if (nlevels(condition) != 2) {
    counts <- paste(names(table(condition)), as.integer(table(condition)), sep = "=", collapse = ", ")
    raw_counts <- paste(names(table(meta$condition_raw, useNA = "ifany")), as.integer(table(meta$condition_raw, useNA = "ifany")), sep = "=", collapse = ", ")
    stop(
      "Contrast ", contrast, " requires two condition levels after metadata curation and filtering. ",
      "Observed binary levels: ", counts, ". Raw condition counts: ", raw_counts, ". ",
      "Inspect results/tables/curated_metadata_all.tsv or run scripts/export_curated_metadata.R.",
      call. = FALSE
    )
  }

  cov_df <- NULL
  if (!is.null(adjust_covariates)) {
    cov_df <- adjust_covariates[meta$sample, , drop = FALSE]
    cov_df <- as.data.frame(scale(cov_df))
    cov_df <- cov_df[, vapply(cov_df, function(x) any(is.finite(x)), logical(1)), drop = FALSE]
    names(cov_df) <- make.names(names(cov_df), unique = TRUE)
    cov_df <- cov_df[meta$sample, , drop = FALSE]
    ok <- stats::complete.cases(cov_df)
    if (!all(ok)) {
      meta <- meta[ok, , drop = FALSE]
      expr <- expr[, meta$sample, drop = FALSE]
      condition <- droplevels(factor(
        ifelse(meta$condition_raw == cf$case_level, cf$case_level, cf$reference_level),
        levels = c(cf$reference_level, cf$case_level)
      ))
      meta$condition <- condition
      meta$condition_binary <- condition
      cov_df <- cov_df[ok, , drop = FALSE]
    }
  }

  patient_available <- "patient_id" %in% names(meta) && sum(!is.na(meta$patient_id) & nzchar(as.character(meta$patient_id))) > 0
  duplicated_patient <- patient_available && any(duplicated(meta$patient_id[!is.na(meta$patient_id) & nzchar(as.character(meta$patient_id))]))
  paired_ok <- isTRUE(paired_design) && duplicated_patient

  if (paired_ok) {
    meta$patient_id <- factor(meta$patient_id)
    design <- stats::model.matrix(~ 0 + condition + patient_id, data = meta)
    colnames(design) <- make.names(colnames(design), unique = TRUE)
    if (!is.null(cov_df) && ncol(cov_df)) design <- cbind(design, cov_df)
    validate_design_expr(design, expr, meta, context = paste0(" for paired ", contrast))
    contrast_mat <- make_paired_contrast_matrix(design, condition, contrast)
    fit <- limma::lmFit(expr, design)
    fit2_unmoderated <- limma::contrasts.fit(fit, contrasts = contrast_mat)
    fit2 <- limma::eBayes(fit2_unmoderated)
    coef_name <- contrast
  } else if (identical(paired_design, "mixed") && duplicated_patient) {
    # Mixed paired/unpaired datasets: use duplicateCorrelation, but never allow
    # missing block IDs. Singleton/unpaired samples get their own block so the
    # block vector remains valid and the design has exactly one row per sample.
    block <- as.character(meta$patient_id)
    block[is.na(block) | !nzchar(block)] <- meta$sample[is.na(block) | !nzchar(block)]
    block <- factor(block)
    design <- stats::model.matrix(~ condition, data = meta)
    colnames(design) <- make.names(colnames(design), unique = TRUE)
    if (!is.null(cov_df) && ncol(cov_df)) design <- cbind(design, cov_df)
    validate_design_expr(design, expr, meta, context = paste0(" for mixed ", contrast))
    coef_name <- coefficient_for_unpaired_model(design, condition, contrast)
    corfit <- limma::duplicateCorrelation(expr, design = design, block = block)
    fit <- limma::lmFit(expr, design = design, block = block, correlation = corfit$consensus.correlation)
    fit2_unmoderated <- fit
    fit2 <- limma::eBayes(fit2_unmoderated)
  } else {
    design <- stats::model.matrix(~ condition, data = meta)
    colnames(design) <- make.names(colnames(design), unique = TRUE)
    if (!is.null(cov_df) && ncol(cov_df)) design <- cbind(design, cov_df)
    validate_design_expr(design, expr, meta, context = paste0(" for unpaired ", contrast))
    coef_name <- coefficient_for_unpaired_model(design, condition, contrast)
    fit <- limma::lmFit(expr, design)
    fit2_unmoderated <- fit
    fit2 <- limma::eBayes(fit2_unmoderated)
  }

  se_vec <- ordinary_se_from_limma(fit2_unmoderated, coef_name)
  tt <- safe_top_table(fit2, coef_name, contrast, design = design) |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      contrast = contrast,
      n_tumor = sum(meta$condition_binary == levels(meta$condition_binary)[2]),
      n_control = sum(meta$condition_binary == levels(meta$condition_binary)[1]),
      se = as.numeric(se_vec[.data$gene]),
      se_source = "ordinary_limma_unmoderated",
      moderated_t = .data$t,
      fdr = stats::p.adjust(.data$P.Value, method = "BH")
    )
  tt
}

run_de_one <- function(ds, cfg, contrast = NULL, adjusted_scores = NULL) {
  contrast <- contrast %||% ds$config$primary_contrast %||% "tumor_vs_control"
  expr_meta <- average_technical_replicates(
    ds$expr,
    ds$meta,
    enabled = isTRUE(ds$config$technical_replicates)
  )

  if (!is.null(adjusted_scores)) {
    adjusted_scores <- as.data.frame(adjusted_scores)
    if (!all(expr_meta$meta$sample %in% rownames(adjusted_scores))) {
      if (isTRUE(ds$config$technical_replicates) && "technical_group" %in% names(ds$meta)) {
        tmp <- adjusted_scores[ds$meta$sample, , drop = FALSE]
        tmp$technical_group <- ds$meta$technical_group
        adjusted_scores <- tmp |>
          dplyr::group_by(.data$technical_group) |>
          dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
          tibble::column_to_rownames("technical_group")
      }
    }
  }

  run_limma_contrast(
    expr = expr_meta$expr,
    meta = expr_meta$meta,
    contrast = contrast,
    paired_design = ds$config$paired_design %||% FALSE,
    adjust_covariates = adjusted_scores
  ) |>
    dplyr::mutate(dataset = ds$dataset, platform = ds$platform, .before = 1)
}

run_all_de <- function(dataset_list, cfg, contrast = "tumor_vs_control") {
  out <- purrr::map(dataset_list, ~ run_de_one(.x, cfg, contrast = contrast))
  names(out) <- names(dataset_list)
  all <- dplyr::bind_rows(out)
  write_tsv(all, file.path(cfg$project$output_dir, "tables", paste0("within_study_de_", contrast, ".tsv")))
  out
}

run_score_de_one <- function(score_matrix, meta, ds_cfg, contrast = "tumor_vs_control") {
  # score_matrix: features x samples
  tmp <- list(expr = score_matrix, meta = meta, config = ds_cfg, dataset = unique(meta$dataset)[1], platform = "score")
  run_de_one(tmp, cfg = list(), contrast = contrast)
}
