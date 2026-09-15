# DepMap and Open Targets target-prioritization layer.
# DepMap is handled from local files to avoid changing portal URLs. Open Targets
# GraphQL is optional and fails gracefully when offline or schema changes.

get_target_cfg <- function(cfg) get_validation_subcfg(cfg, "targets")

candidate_genes_for_targeting <- function(cfg) {
  tabs <- load_atlas_signature_tables(cfg)
  hubs <- tabs$main_hubs
  if (!nrow(hubs)) hubs <- tabs$strict_hubs
  gcol <- intersect(c("feature", "gene"), colnames(hubs))[1]
  if (is.na(gcol)) return(character())
  n <- get_target_cfg(cfg)$top_n_candidates %||% 100
  unique(utils::head(normalise_gene_symbols(hubs[[gcol]]), n))
}

depmap_column_symbol <- function(cols) {
  # DepMap gene columns commonly look like "KRAS (3845)" or
  # "KRAS (ENSG00000133703)". Keep the left-side symbol for matching while
  # leaving non-gene ID columns unmatched.
  sym <- safe_to_utf8(cols)
  sym <- sub("\\s*\\(.*$", "", sym)
  sym <- sub("\\..*$", "", sym)
  normalise_gene_symbols(sym)
}

find_depmap_gene_columns <- function(tbl, genes) {
  if (is.null(tbl) || !length(genes)) return(stats::setNames(character(), character()))
  cols <- colnames(tbl)
  gene_cols <- setdiff(cols, c("depmap_id", "DepMap_ID", "ModelID", "model_id", "ProfileID", "profile_id", "cell_line", "CellLineName"))
  symbols <- depmap_column_symbol(gene_cols)
  genes_norm <- unique(normalise_gene_symbols(genes))
  map <- purrr::map_chr(genes_norm, function(g) {
    idx <- which(symbols == g)
    if (length(idx)) return(gene_cols[idx[1]])
    NA_character_
  })
  names(map) <- genes_norm
  map[!is.na(map)]
}

depmap_sample_ids_from_info <- function(sample_info) {
  if (is.null(sample_info) || !nrow(sample_info)) return(character())
  id_col <- intersect(c("DepMap_ID", "depmap_id", "ModelID", "model_id"), colnames(sample_info))[1]
  if (is.na(id_col)) id_col <- colnames(sample_info)[1]
  unique(safe_to_utf8(sample_info[[id_col]]))
}

find_existing_depmap_file <- function(primary, patterns = character()) {
  if (!is.null(primary) && nzchar(primary) && file.exists(primary)) return(primary)
  dirs <- unique(dirname(c(primary %||% "", patterns)))
  dirs <- dirs[!is.na(dirs) & nzchar(dirs) & dir.exists(dirs)]
  hits <- character()
  for (d in dirs) {
    files <- list.files(d, full.names = TRUE, recursive = FALSE)
    for (pat in basename(patterns)) {
      hits <- c(hits, files[grepl(pat, basename(files), ignore.case = TRUE)])
    }
  }
  hits <- unique(hits[file.exists(hits)])
  if (length(hits)) hits[1] else primary
}

read_depmap_matrix <- function(path, sample_info = NULL) {
  x <- read_table_auto(path)
  ids <- depmap_sample_ids_from_info(sample_info)
  id_candidates <- c("DepMap_ID", "depmap_id", "ModelID", "model_id", "ProfileID", "profile_id", "cell_line", "CellLineName")

  # Standard model-by-gene layout: one row per model, gene columns such as
  # "KRAS (3845)" or "KRAS (ENSG00000133703)".
  id_col <- intersect(id_candidates, colnames(x))[1]
  if (!is.na(id_col)) {
    return(x |> dplyr::rename(depmap_id = dplyr::all_of(id_col)))
  }

  # Some DepMap omics exports are gene-by-model: first column is the gene
  # feature and the remaining columns are DepMap model IDs. Transpose these.
  model_cols <- if (length(ids)) intersect(colnames(x), ids) else character()
  if (length(model_cols) >= 5) {
    gene_col <- colnames(x)[1]
    genes <- safe_to_utf8(x[[gene_col]])
    keep <- !is.na(genes) & nzchar(genes)
    if (sum(keep) >= 5) {
      mat <- as.data.frame(t(as.matrix(x[keep, model_cols, drop = FALSE])), check.names = FALSE)
      colnames(mat) <- make.unique(genes[keep])
      mat <- tibble::rownames_to_column(mat, "depmap_id")
      return(clean_character_columns(mat))
    }
  }

  # Fallback: first column may still be model IDs but named differently.
  first_col <- colnames(x)[1]
  first_vals <- safe_to_utf8(x[[first_col]])
  if (!length(ids) || sum(first_vals %in% ids, na.rm = TRUE) >= 5 || grepl("ACH-", paste(utils::head(first_vals, 20), collapse = " "))) {
    return(x |> dplyr::rename(depmap_id = dplyr::all_of(first_col)))
  }

  # Last resort: keep a depmap_id column so downstream joins fail gracefully and
  # the diagnostics report "not_found" instead of crashing.
  x |> dplyr::rename(depmap_id = dplyr::all_of(first_col))
}

load_depmap_inputs <- function(cfg) {
  tcfg <- get_target_cfg(cfg)
  sample_info <- if (file.exists(tcfg$depmap_sample_info_file %||% "")) read_table_auto(tcfg$depmap_sample_info_file) else NULL
  expr_file <- find_existing_depmap_file(
    tcfg$depmap_expression_file %||% "",
    patterns = file.path(dirname(tcfg$depmap_expression_file %||% "data/external/depmap/x"), c(
      "OmicsExpression.*\\.csv(\\.gz)?$",
      "Expression.*TPM.*\\.csv(\\.gz)?$",
      ".*Expression.*ProteinCoding.*\\.csv(\\.gz)?$"
    ))
  )
  crispr_file <- find_existing_depmap_file(
    tcfg$depmap_crispr_gene_effect_file %||% "",
    patterns = file.path(dirname(tcfg$depmap_crispr_gene_effect_file %||% "data/external/depmap/x"), c(
      "CRISPR.*Gene.*Effect.*\\.csv(\\.gz)?$",
      ".*GeneEffect.*\\.csv(\\.gz)?$"
    ))
  )
  list(
    expression = if (!is.null(expr_file) && nzchar(expr_file) && file.exists(expr_file)) read_depmap_matrix(expr_file, sample_info = sample_info) else NULL,
    crispr = if (!is.null(crispr_file) && nzchar(crispr_file) && file.exists(crispr_file)) read_depmap_matrix(crispr_file, sample_info = sample_info) else NULL,
    sample_info = sample_info,
    expression_file = expr_file,
    crispr_file = crispr_file
  )
}

filter_pancreas_depmap <- function(sample_info) {
  if (is.null(sample_info) || !nrow(sample_info)) return(NULL)
  sample_info <- clean_character_columns(sample_info)
  names(sample_info) <- make.names(names(sample_info), unique = TRUE)
  id_col <- intersect(c("DepMap_ID", "depmap_id", "ModelID", "model_id"), names(sample_info))[1]
  if (is.na(id_col)) id_col <- names(sample_info)[1]
  sample_info$depmap_id <- safe_to_utf8(sample_info[[id_col]])
  collapsed <- collapse_rows_lower(sample_info)
  sample_info$lineage_pancreas_flag <- grepl("pancre", collapsed)
  sample_info |> dplyr::filter(.data$lineage_pancreas_flag)
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

summarise_depmap_gene_set <- function(depmap, genes, cfg = NULL) {
  # Always return a table with a gene column, even when no DepMap gene columns
  # match the candidate symbols. This avoids dplyr join failures and makes the
  # report explicit about missing DepMap coverage.
  genes <- unique(normalise_gene_symbols(genes))
  base <- tibble::tibble(gene = genes)
  if (!length(genes)) return(base)

  sample_info <- filter_pancreas_depmap(depmap$sample_info)
  pancreas_ids <- if (!is.null(sample_info)) sample_info$depmap_id else character()

  expr_res <- tibble::tibble(gene = character())
  if (!is.null(depmap$expression)) {
    map <- find_depmap_gene_columns(depmap$expression, genes)
    expr_res <- purrr::map_dfr(names(map), function(g) {
      col <- map[[g]]
      vals_all <- suppressWarnings(as.numeric(depmap$expression[[col]]))
      vals_pan <- vals_all[depmap$expression$depmap_id %in% pancreas_ids]
      mean_pan <- safe_mean(vals_pan)
      tibble::tibble(
        gene = g,
        depmap_expression_col = col,
        n_all_expression = sum(is.finite(vals_all)),
        n_pancreas_expression = sum(is.finite(vals_pan)),
        mean_expression_all = safe_mean(vals_all),
        mean_expression_pancreas = mean_pan,
        pancreas_expression_percentile = if (is.finite(mean_pan)) mean(vals_all <= mean_pan, na.rm = TRUE) else NA_real_
      )
    })
  }

  crispr_res <- tibble::tibble(gene = character())
  if (!is.null(depmap$crispr)) {
    map <- find_depmap_gene_columns(depmap$crispr, genes)
    crispr_res <- purrr::map_dfr(names(map), function(g) {
      col <- map[[g]]
      vals_all <- suppressWarnings(as.numeric(depmap$crispr[[col]]))
      vals_pan <- vals_all[depmap$crispr$depmap_id %in% pancreas_ids]
      tibble::tibble(
        gene = g,
        depmap_crispr_col = col,
        n_all_crispr = sum(is.finite(vals_all)),
        n_pancreas_crispr = sum(is.finite(vals_pan)),
        median_gene_effect_all = safe_median(vals_all),
        median_gene_effect_pancreas = safe_median(vals_pan),
        n_pancreas_strong_dependency_lt_minus1 = sum(vals_pan < -1, na.rm = TRUE),
        n_pancreas_dependency_lt_minus05 = sum(vals_pan < -0.5, na.rm = TRUE)
      )
    })
  }

  dep <- base |>
    dplyr::left_join(expr_res, by = "gene") |>
    dplyr::left_join(crispr_res, by = "gene")

  needed_depmap_cols <- c(
    "depmap_expression_col", "n_all_expression", "n_pancreas_expression",
    "mean_expression_all", "mean_expression_pancreas", "pancreas_expression_percentile",
    "depmap_crispr_col", "n_all_crispr", "n_pancreas_crispr",
    "median_gene_effect_all", "median_gene_effect_pancreas",
    "n_pancreas_strong_dependency_lt_minus1", "n_pancreas_dependency_lt_minus05"
  )
  for (nm in needed_depmap_cols) {
    if (!nm %in% colnames(dep)) {
      dep[[nm]] <- if (grepl("_col$", nm)) NA_character_ else NA_real_
    }
  }

  dep <- dep |>
    dplyr::mutate(
      depmap_expression_status = dplyr::if_else(!is.na(.data$depmap_expression_col), "matched", "not_found"),
      depmap_crispr_status = dplyr::if_else(!is.na(.data$depmap_crispr_col), "matched", "not_found")
    )
  if (!is.null(cfg)) {
    write_tsv_safe(dep, file.path(external_table_dir(cfg), "depmap_candidate_summary.tsv"))
  }
  dep
}

symbol_to_ensembl <- function(genes) {
  genes <- normalise_gene_symbols(genes)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(tibble::tibble(gene = genes, ensembl = NA_character_))
  }
  ann <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("SYMBOL", "ENSEMBL")) |>
    dplyr::transmute(gene = normalise_gene_symbols(.data$SYMBOL), ensembl = .data$ENSEMBL) |>
    dplyr::filter(!is.na(.data$ensembl)) |>
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  tibble::tibble(gene = genes) |> dplyr::left_join(ann, by = "gene")
}

opentargets_query_target <- function(ensembl_id, disease_id = NULL, timeout_sec = 30) {
  if (!requireNamespace("httr2", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  # v3.11: Open Targets Platform API 26.x no longer exposes target.knownDrugs.
  # Use drugAndClinicalCandidates and keep associatedDiseases + tractability.
  query <- '
  query targetInfo($ensemblId: String!) {
    target(ensemblId: $ensemblId) {
      id
      approvedSymbol
      approvedName
      biotype
      tractability {
        modality
        label
        value
      }
      associatedDiseases(page: {index: 0, size: 50}, orderByScore: "score desc") {
        count
        rows {
          score
          disease { id name }
        }
      }
      drugAndClinicalCandidates {
        count
        rows {
          maxClinicalStage
          drug { id name }
          diseases {
            diseaseFromSource
            disease { id name }
          }
        }
      }
    }
  }'
  req <- httr2::request("https://api.platform.opentargets.org/api/v4/graphql") |>
    httr2::req_method("POST") |>
    httr2::req_timeout(timeout_sec) |>
    httr2::req_body_json(list(query = query, variables = list(ensemblId = ensembl_id))) |>
    httr2::req_error(is_error = function(resp) FALSE)
  resp <- tryCatch(httr2::req_perform(req), error = function(e) e)
  if (inherits(resp, "error")) return(list(error = conditionMessage(resp)))
  status <- httr2::resp_status(resp)
  dat <- tryCatch(httr2::resp_body_json(resp), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(dat$errors)) {
    err <- paste(vapply(dat$errors, function(x) x$message %||% as.character(x), character(1)), collapse = " | ")
    return(list(error = paste0("GraphQL error: ", err)))
  }
  if (status < 200 || status >= 300) {
    return(list(error = paste0("HTTP ", status, " ", httr2::resp_status_desc(resp))))
  }
  dat
}

summarise_opentargets_response <- function(gene, ensembl, dat, disease_id = NULL, disease_regex = "pancre") {
  if (is.null(dat) || !is.null(dat$error)) {
    return(tibble::tibble(gene = gene, ensembl = ensembl, opentargets_status = "error_or_skipped", opentargets_error = dat$error %||% NA_character_))
  }
  target <- dat$data$target
  if (is.null(target)) return(tibble::tibble(gene = gene, ensembl = ensembl, opentargets_status = "no_target"))
  tr <- target$tractability %||% list()
  tract <- if (length(tr)) paste(unique(purrr::map_chr(tr, function(x) paste(x$modality %||% NA, x$label %||% NA, x$value %||% NA, sep = ":"))), collapse = "; ") else NA_character_
  dis <- target$associatedDiseases$rows %||% list()
  dis_tbl <- if (length(dis)) purrr::map_dfr(dis, function(x) tibble::tibble(disease_id = x$disease$id %||% NA_character_, disease_name = x$disease$name %||% NA_character_, disease_score = x$score %||% NA_real_)) else tibble::tibble()
  has_disease_id <- !is.null(disease_id) && length(disease_id) == 1 && !is.na(disease_id) && nzchar(as.character(disease_id))
  if (nrow(dis_tbl)) {
    if (has_disease_id) {
      pdac <- dis_tbl |>
        dplyr::filter(.data$disease_id == !!as.character(disease_id) | grepl(disease_regex, .data$disease_name, ignore.case = TRUE)) |>
        dplyr::arrange(dplyr::desc(.data$disease_score))
    } else {
      pdac <- dis_tbl |>
        dplyr::filter(grepl(disease_regex, .data$disease_name, ignore.case = TRUE)) |>
        dplyr::arrange(dplyr::desc(.data$disease_score))
    }
    top_dis <- dis_tbl |> dplyr::arrange(dplyr::desc(.data$disease_score)) |> utils::head(1)
  } else {
    pdac <- tibble::tibble(); top_dis <- tibble::tibble()
  }
  kd <- target$drugAndClinicalCandidates$rows %||% target$knownDrugs$rows %||% list()
  known <- if (length(kd)) {
    paste(unique(purrr::map_chr(kd, function(x) {
      drug_name <- x$drug$name %||% NA_character_
      stage <- x$maxClinicalStage %||% x$phase %||% NA_character_
      dis <- x$diseases %||% list()
      dis_name <- if (length(dis)) {
        paste(unique(purrr::map_chr(dis, function(d) d$disease$name %||% d$diseaseFromSource %||% NA_character_)), collapse = ",")
      } else {
        x$disease$name %||% NA_character_
      }
      paste(drug_name, dis_name, stage, sep = ":")
    })), collapse = "; ")
  } else NA_character_
  pancreas_known_vec <- if (length(kd)) {
    purrr::map_chr(kd, function(x) {
      drug_name <- x$drug$name %||% NA_character_
      stage <- x$maxClinicalStage %||% x$phase %||% NA_character_
      dis <- x$diseases %||% list()
      dis_name <- if (length(dis)) {
        paste(unique(purrr::map_chr(dis, function(d) d$disease$name %||% d$diseaseFromSource %||% NA_character_)), collapse = ",")
      } else {
        x$disease$name %||% NA_character_
      }
      if (grepl(disease_regex, dis_name, ignore.case = TRUE)) paste(drug_name, dis_name, stage, sep = ":") else NA_character_
    })
  } else character()
  pancreas_known_vec <- unique(pancreas_known_vec[!is.na(pancreas_known_vec) & nzchar(pancreas_known_vec)])
  pancreas_known <- if (length(pancreas_known_vec)) paste(pancreas_known_vec, collapse = "; ") else NA_character_
  known_count <- target$drugAndClinicalCandidates$count %||% target$knownDrugs$count %||% length(kd)
  tibble::tibble(
    gene = gene,
    ensembl = ensembl,
    opentargets_status = "ok",
    approved_symbol = target$approvedSymbol %||% NA_character_,
    approved_name = target$approvedName %||% NA_character_,
    biotype = target$biotype %||% NA_character_,
    tractability = tract,
    top_disease_name = if (nrow(top_dis)) top_dis$disease_name[1] else NA_character_,
    top_disease_score = if (nrow(top_dis)) top_dis$disease_score[1] else NA_real_,
    pancreas_disease_name = if (nrow(pdac)) pdac$disease_name[1] else NA_character_,
    pancreas_disease_score = if (nrow(pdac)) pdac$disease_score[1] else NA_real_,
    known_drug_count = suppressWarnings(as.numeric(known_count)),
    known_drugs = known,
    pancreas_known_drug_count = length(pancreas_known_vec),
    pancreas_known_drugs = pancreas_known
  )
}

query_opentargets_for_candidates <- function(cfg, genes) {
  tcfg <- get_target_cfg(cfg)
  if (!isTRUE(tcfg$query_opentargets %||% FALSE)) {
    res <- symbol_to_ensembl(genes) |> dplyr::mutate(opentargets_status = "skipped_config_query_opentargets_false")
    write_tsv_safe(res, file.path(external_table_dir(cfg), "opentargets_annotation.tsv"))
    return(res)
  }
  ann <- symbol_to_ensembl(genes)
  res <- purrr::pmap_dfr(ann, function(gene, ensembl) {
    if (is.na(ensembl)) return(tibble::tibble(gene = gene, ensembl = NA_character_, opentargets_status = "no_ensembl"))
    Sys.sleep(tcfg$opentargets_sleep_sec %||% 0.2)
    dat <- opentargets_query_target(ensembl, disease_id = tcfg$opentargets_disease_id %||% NULL, timeout_sec = tcfg$opentargets_timeout_sec %||% 30)
    summarise_opentargets_response(gene, ensembl, dat, disease_id = tcfg$opentargets_disease_id %||% NULL)
  })
  write_tsv_safe(res, file.path(external_table_dir(cfg), "opentargets_annotation.tsv"))
  res
}


# v3.13: after joining strict-hub, manual, DepMap, and Open Targets tables,
# several annotation columns can exist as .x/.y suffixed duplicates. Coalesce
# them back into canonical columns before computing targetability. This prevents
# manually curated tractability annotations (e.g. SLC6A14 transporter / plasma
# membrane, LY75 receptor) from being silently ignored.
coalesce_target_column_variants <- function(tbl, base_name) {
  variants <- intersect(
    c(base_name, paste0(base_name, ".x"), paste0(base_name, ".y"),
      paste0(base_name, ".x.x"), paste0(base_name, ".x.y"),
      paste0(base_name, ".y.x"), paste0(base_name, ".y.y")),
    colnames(tbl)
  )
  if (!length(variants)) return(tbl)
  vals <- tbl[[variants[1]]]
  if (length(variants) > 1) {
    for (v in variants[-1]) {
      vals <- dplyr::coalesce(vals, tbl[[v]])
    }
  }
  tbl[[base_name]] <- vals
  tbl
}

collapse_target_annotation_columns <- function(tbl) {
  for (nm in c(
    "target_class", "tractability_manual", "tractability",
    "subcellular_location", "annotation_source", "annotation_bonus"
  )) {
    tbl <- coalesce_target_column_variants(tbl, nm)
  }
  tbl
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(is.finite(x) & x != 0)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "yes", "y", "1")
}


make_target_prioritization_table <- function(cfg) {
  genes <- candidate_genes_for_targeting(cfg)
  depmap <- load_depmap_inputs(cfg)
  dep <- if (!is.null(depmap$expression) || !is.null(depmap$crispr)) summarise_depmap_gene_set(depmap, genes, cfg = cfg) else tibble::tibble(gene = genes)
  ot <- query_opentargets_for_candidates(cfg, genes)
  strict <- load_atlas_signature_tables(cfg)$strict_hubs
  if (nrow(strict)) {
    gcol <- intersect(c("feature", "gene"), colnames(strict))[1]
    if (!is.na(gcol)) {
      strict <- strict |> dplyr::mutate(gene = normalise_gene_symbols(.data[[gcol]]))
    }
  }
  if (!"gene" %in% colnames(strict)) strict <- tibble::tibble(gene = character())

  drug_manual <- file.path("resources", "druggability_manual_seed.csv")
  manual <- if (file.exists(drug_manual)) {
    readr::read_csv(drug_manual, show_col_types = FALSE) |>
      dplyr::mutate(gene = normalise_gene_symbols(.data$gene))
  } else {
    tibble::tibble(gene = genes)
  }
  if (!"gene" %in% colnames(manual)) manual <- tibble::tibble(gene = character())
  if (!"gene" %in% colnames(dep)) dep <- tibble::tibble(gene = genes)
  if (!"gene" %in% colnames(ot)) ot <- tibble::tibble(gene = genes)

  pri <- tibble::tibble(gene = genes) |>
    dplyr::left_join(strict, by = "gene") |>
    dplyr::left_join(manual, by = "gene") |>
    dplyr::left_join(dep, by = "gene") |>
    dplyr::left_join(ot, by = "gene") |>
    collapse_target_annotation_columns()
  needed_numeric <- c(
    "strict_hub_score", "median_gene_effect_pancreas", "mean_expression_pancreas",
    "mean_expression_all", "pancreas_disease_score", "pancreas_known_drug_count",
    "n_pancreas_crispr", "n_pancreas_dependency_lt_minus05", "annotation_bonus"
  )
  needed_character <- c("target_class", "tractability_manual", "tractability", "subcellular_location", "tme_class")
  for (nm in needed_numeric) if (!nm %in% colnames(pri)) pri[[nm]] <- NA_real_
  for (nm in needed_character) if (!nm %in% colnames(pri)) pri[[nm]] <- NA_character_
  if (!"retained_after_tme_adjustment" %in% colnames(pri)) pri$retained_after_tme_adjustment <- FALSE
  pri <- pri |>
    dplyr::mutate(
      retained_for_tumor_associated_program = as_logical_safe(.data$retained_after_tme_adjustment %||% FALSE) |
        (.data$tme_class %in% c("tumor_associated_retained_under_S1", "tumor_intrinsic_retained")),
      target_annotation_text = paste(
        .data$target_class, .data$tractability_manual,
        .data$tractability, .data$subcellular_location,
        sep = " "
      ),
      targetability_bonus = dplyr::case_when(
        grepl("small_molecule|small molecule|transporter|antibody|ligand|lectin|receptor", .data$target_annotation_text, ignore.case = TRUE) ~ 2,
        grepl("plasma membrane|extracellular|Approved Drug:TRUE|Advanced Clinical:TRUE|Phase 1 Clinical:TRUE|High-Quality Ligand:TRUE|Druggable Family:TRUE", .data$target_annotation_text, ignore.case = TRUE) ~ 1,
        TRUE ~ 0
      ),
      pancreas_dependency_fraction_lt_minus05 = dplyr::if_else(
        is.finite(.data$n_pancreas_crispr) & .data$n_pancreas_crispr > 0,
        .data$n_pancreas_dependency_lt_minus05 / .data$n_pancreas_crispr,
        NA_real_
      ),
      common_essential_dependency_flag = is.finite(.data$median_gene_effect_pancreas) &
        .data$median_gene_effect_pancreas < -0.5 &
        is.finite(.data$pancreas_dependency_fraction_lt_minus05) &
        .data$pancreas_dependency_fraction_lt_minus05 >= 0.50,
      depmap_dependency_bonus = dplyr::case_when(
        is.finite(.data$median_gene_effect_pancreas) & .data$median_gene_effect_pancreas < -0.5 ~ 2,
        is.finite(.data$median_gene_effect_pancreas) & .data$median_gene_effect_pancreas < -0.25 ~ 1,
        TRUE ~ 0
      ),
      expression_bonus = dplyr::case_when(
        is.finite(.data$mean_expression_pancreas) & is.finite(.data$mean_expression_all) & .data$mean_expression_pancreas > .data$mean_expression_all ~ 1,
        TRUE ~ 0
      ),
      opentargets_bonus = dplyr::case_when(
        is.finite(.data$pancreas_disease_score) & .data$pancreas_disease_score > 0.2 ~ 2,
        is.finite(.data$pancreas_disease_score) & .data$pancreas_disease_score > 0 ~ 1,
        is.finite(.data$pancreas_known_drug_count) & .data$pancreas_known_drug_count > 0 ~ 1,
        TRUE ~ 0
      ),
      external_target_priority_score = dplyr::coalesce(.data$strict_hub_score, 0) +
        .data$targetability_bonus + .data$depmap_dependency_bonus +
        .data$expression_bonus + .data$opentargets_bonus,
      candidate_lane = dplyr::case_when(
        .data$retained_for_tumor_associated_program & .data$targetability_bonus > 0 ~ "retained_tumor_associated_tractable",
        .data$retained_for_tumor_associated_program ~ "retained_tumor_associated",
        .data$common_essential_dependency_flag ~ "common_essential_dependency",
        .data$tme_class %in% c("composition_covarying", "microenvironment_sensitive") ~ "composition_covarying_marker",
        TRUE ~ "exploratory_or_low_priority"
      ),
      pdac_specific_annotation_bonus = dplyr::case_when(
        is.finite(.data$pancreas_disease_score) & .data$pancreas_disease_score > 0.2 ~ 2,
        is.finite(.data$pancreas_disease_score) & .data$pancreas_disease_score > 0 ~ 1,
        TRUE ~ 0
      ),
      main_text_candidate_score = dplyr::coalesce(.data$strict_hub_score, 0) +
        3 * as.numeric(.data$retained_for_tumor_associated_program) +
        .data$targetability_bonus + .data$expression_bonus +
        .data$pdac_specific_annotation_bonus -
        2 * as.numeric(.data$common_essential_dependency_flag) -
        1 * as.numeric(.data$tme_class %in% c("composition_covarying", "microenvironment_sensitive"))
    ) |>
    dplyr::arrange(dplyr::desc(.data$main_text_candidate_score), dplyr::desc(.data$external_target_priority_score))
  write_tsv_safe(pri, file.path(external_table_dir(cfg), "external_target_prioritization.tsv"))
  if (nrow(pri)) {
    p <- pri |> utils::head(30) |>
      ggplot2::ggplot(ggplot2::aes(
        x = stats::reorder(gene, main_text_candidate_score),
        y = main_text_candidate_score,
        fill = candidate_lane
      )) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Gene",
        y = "Conservative main-text candidate score",
        fill = "Candidate lane",
        title = "Conservative DepMap/Open Targets target annotation"
      ) +
      ggplot2::theme_bw(base_size = 9)
    ggplot2::ggsave(file.path(external_figure_dir(cfg), "external_target_prioritization.png"), p, width = 8, height = 8, dpi = 300)
  }
  pri
}
