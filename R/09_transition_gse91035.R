
# Defensive fallback for standalone sourcing. The canonical definition lives in R/00_utils.R.
if (!exists("normalise_gene_symbols", mode = "function")) {
  normalise_gene_symbols <- function(x) {
    x <- as.character(x)
    x <- sub("^.*\\|", "", x)
    x <- sub("\\..*$", "", x)
    toupper(trimws(x))
  }
}

# Normal-benign-PDAC transition analysis -------------------------------------

run_gse91035_transition <- function(dataset_list, cfg) {
  if (!"GSE91035" %in% names(dataset_list)) return(tibble::tibble())
  ds <- dataset_list[["GSE91035"]]
  meta <- ds$meta |> dplyr::filter(.data$condition %in% c("control", "benign", "tumor"))
  if (length(unique(meta$condition)) < 3) {
    warning("GSE91035 does not have all three classes after curation; transition analysis skipped.")
    return(tibble::tibble())
  }
  expr <- ds$expr[, meta$sample, drop = FALSE]
  meta <- meta |>
    dplyr::mutate(
      transition_stage = dplyr::case_when(
        .data$condition == "control" ~ 0,
        .data$condition == "benign" ~ 1,
        .data$condition == "tumor" ~ 2,
        TRUE ~ NA_real_
      ),
      condition3 = factor(.data$condition, levels = c("control", "benign", "tumor"))
    )

  design_trend <- model.matrix(~ transition_stage, data = meta)
  fit_trend <- limma::eBayes(limma::lmFit(expr, design_trend))
  trend <- limma::topTable(fit_trend, coef = "transition_stage", number = Inf, sort.by = "none") |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(fdr_trend = bh(.data$P.Value), trend_logFC_per_stage = .data$logFC)

  design_anova <- stats::model.matrix(~ 0 + condition3, data = meta)
  fit_anova <- limma::eBayes(limma::lmFit(expr, design_anova))
  anova_res <- limma::topTable(fit_anova, coef = seq_len(ncol(design_anova)), number = Inf, sort.by = "none") |>
    tibble::rownames_to_column("gene") |>
    dplyr::transmute(gene = .data$gene, p_anova = .data$P.Value, fdr_anova = bh(.data$P.Value), F_anova = .data$F)

  pairwise <- list(
    pdac_vs_normal = run_limma_contrast(expr, meta |> dplyr::mutate(condition = dplyr::recode(.data$condition, control = "control", tumor = "tumor", benign = "benign")), "pdac_vs_normal", paired_design = FALSE),
    pdac_vs_benign = run_limma_contrast(expr, meta, "pdac_vs_benign", paired_design = FALSE),
    benign_vs_normal = run_limma_contrast(expr, meta, "benign_vs_normal", paired_design = FALSE)
  )

  out <- trend |>
    dplyr::select(.data$gene, .data$trend_logFC_per_stage, .data$P.Value, .data$fdr_trend) |>
    dplyr::rename(p_trend = .data$P.Value) |>
    dplyr::left_join(pairwise$pdac_vs_benign |> dplyr::select(.data$gene, pdac_vs_benign_logFC = .data$logFC, pdac_vs_benign_fdr = .data$fdr), by = "gene") |>
    dplyr::left_join(pairwise$benign_vs_normal |> dplyr::select(.data$gene, benign_vs_normal_logFC = .data$logFC, benign_vs_normal_fdr = .data$fdr), by = "gene") |>
    dplyr::left_join(pairwise$pdac_vs_normal |> dplyr::select(.data$gene, pdac_vs_normal_logFC = .data$logFC, pdac_vs_normal_fdr = .data$fdr), by = "gene") |>
    dplyr::left_join(anova_res, by = "gene") |>
    dplyr::mutate(
      nonmonotonic_screen_class = dplyr::case_when(
        fdr_anova < (cfg$analysis$transition_program$fdr %||% 0.05) &
          !(fdr_trend < (cfg$analysis$transition_program$fdr %||% 0.05)) ~ "nonmonotonic_candidate",
        fdr_anova < (cfg$analysis$transition_program$fdr %||% 0.05) ~ "three_class_difference",
        TRUE ~ "not_three_class_significant"
      ),
      monotonic_class = dplyr::case_when(
        fdr_trend < (cfg$analysis$transition_program$fdr %||% 0.05) &
          pdac_vs_benign_fdr < 0.05 & abs(pdac_vs_benign_logFC) >= (cfg$analysis$transition_program$abs_logfc_pdac_vs_benign %||% 0.4) ~ "benign_to_pdac_transition",
        benign_vs_normal_fdr < 0.05 & sign(benign_vs_normal_logFC) == sign(pdac_vs_normal_logFC) ~ "early_benign_shift",
        TRUE ~ "not_transition_core"
      )
    ) |>
    dplyr::arrange(.data$fdr_trend, dplyr::desc(abs(.data$trend_logFC_per_stage)))

  write_tsv(out, file.path(cfg$project$output_dir, "tables", "gse91035_normal_benign_pdac_transition.tsv"))
  write_tsv(dplyr::bind_rows(pairwise, .id = "contrast"), file.path(cfg$project$output_dir, "tables", "gse91035_pairwise_de.tsv"))
  out
}


annotate_transition_identifier_class <- function(tbl, cfg = NULL) {
  if (!nrow(tbl)) return(tbl)
  gene_col <- intersect(c("gene", "feature", "symbol"), colnames(tbl))[1]
  if (is.na(gene_col)) return(tbl)
  sy <- normalise_gene_symbols(tbl[[gene_col]])
  mapped_entrez <- rep(NA_character_, length(sy))
  if (quiet_require("AnnotationDbi") && quiet_require("org.Hs.eg.db")) {
    mapped_entrez <- tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = sy, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"),
      error = function(e) rep(NA_character_, length(sy))
    )
    mapped_entrez <- as.character(mapped_entrez[sy])
  }
  noncoding_like <- grepl("^(LOC|LINC|MIR|SNOR|SNHG|RP[0-9]|AL[0-9]|AC[0-9])", sy) |
    grepl("(-AS[0-9]*$|P[0-9]*$|PSEUDO|NCRNA)", sy)
  tbl |>
    dplyr::mutate(
      feature = .data[[gene_col]],
      feature_symbol_normalised = sy,
      hgnc_symbol_mapped_proxy = !is.na(mapped_entrez) & nzchar(mapped_entrez),
      entrez_id_proxy = mapped_entrez,
      noncoding_or_legacy_symbol_like = noncoding_like,
      transition_external_validation_tier = dplyr::case_when(
        .data$hgnc_symbol_mapped_proxy & !.data$noncoding_or_legacy_symbol_like ~ "coding_symbol_proxy",
        .data$noncoding_or_legacy_symbol_like ~ "noncoding_or_legacy_symbol_like",
        TRUE ~ "unresolved_symbol"
      )
    )
}

write_transition_identifier_audit <- function(prog, cfg) {
  if (!nrow(prog)) return(invisible(tibble::tibble()))
  audited <- annotate_transition_identifier_class(prog, cfg)
  summary <- audited |>
    dplyr::count(.data$transition_external_validation_tier, name = "n_genes") |>
    dplyr::mutate(fraction = .data$n_genes / sum(.data$n_genes))
  write_tsv(audited, file.path(cfg$project$output_dir, "tables", "transition_program_identifier_audit.tsv"))
  write_tsv(summary, file.path(cfg$project$output_dir, "tables", "transition_program_identifier_summary.tsv"))
  coding <- audited |>
    dplyr::filter(.data$transition_external_validation_tier == "coding_symbol_proxy")
  write_tsv(coding, file.path(cfg$project$output_dir, "tables", "program_benign_to_pdac_transition_coding_symbol_proxy.tsv"))
  invisible(summary)
}

define_transition_program <- function(transition_res, cfg) {
  if (!nrow(transition_res)) return(tibble::tibble())
  prog <- transition_res |>
    dplyr::filter(.data$monotonic_class == "benign_to_pdac_transition") |>
    dplyr::mutate(
      feature = .data$gene,
      program = dplyr::if_else(.data$trend_logFC_per_stage > 0, "transition_up", "transition_down"),
      transition_model = "ordinal_linear_trend_control_benign_tumor",
      manuscript_placement = "supplementary_until_external_gene_coverage_ge_70pct"
    ) |>
    dplyr::arrange(.data$fdr_trend)
  write_transition_identifier_audit(prog, cfg)
  prog
}
