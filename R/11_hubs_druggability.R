# Integrative hub-gene and druggability annotation ----------------------------

map_symbol_to_ensembl <- function(symbols) {
  symbols <- unique(toupper(symbols))
  if (!quiet_require("org.Hs.eg.db")) return(tibble::tibble(symbol = symbols, ensembl_gene_id = NA_character_))
  ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
  tibble::tibble(symbol = names(ids), ensembl_gene_id = as.character(ids))
}

opentargets_query_one <- function(ensembl_id) {
  endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
  qry <- '
  query target($ensemblId: String!) {
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
      knownDrugs(size: 10) {
        count
        rows {
          phase
          status
          drug { name maximumClinicalTrialPhase }
          disease { name }
        }
      }
    }
  }'
  resp <- httr2::request(endpoint) |>
    httr2::req_body_json(list(query = qry, variables = list(ensemblId = ensembl_id))) |>
    httr2::req_perform()
  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

annotate_opentargets <- function(symbols, cfg) {
  cache_path <- file.path(cfg$project$output_dir, "objects", "opentargets_annotation_cache.rds")
  if (file.exists(cache_path)) return(readRDS(cache_path))
  if (!isTRUE(cfg$analysis$optional_online_annotations)) {
    return(map_symbol_to_ensembl(symbols) |> dplyr::mutate(open_targets_status = "not_queried_optional_online_annotations_false"))
  }
  ids <- map_symbol_to_ensembl(symbols) |> dplyr::filter(!is.na(.data$ensembl_gene_id))
  out <- purrr::map_dfr(seq_len(nrow(ids)), function(i) {
    sym <- ids$symbol[[i]]; ens <- ids$ensembl_gene_id[[i]]
    dat <- tryCatch(opentargets_query_one(ens), error = function(e) NULL)
    if (is.null(dat) || is.null(dat$data$target)) {
      return(tibble::tibble(symbol = sym, ensembl_gene_id = ens, open_targets_status = "query_failed"))
    }
    target <- dat$data$target
    tract <- target$tractability
    tract_txt <- if (length(tract)) paste(vapply(tract, function(x) paste(x$modality, x$label, x$value, sep = ":"), character(1)), collapse = " | ") else NA_character_
    drugs <- target$knownDrugs
    known_drug_count <- drugs$count %||% 0
    top_drugs <- if (!is.null(drugs$rows) && length(drugs$rows)) {
      paste(vapply(drugs$rows, function(x) paste0(x$drug$name %||% NA_character_, " (phase ", x$phase %||% NA_character_, "; ", x$disease$name %||% NA_character_, ")"), character(1)), collapse = " | ")
    } else NA_character_
    tibble::tibble(
      symbol = sym,
      ensembl_gene_id = ens,
      approved_name = target$approvedName %||% NA_character_,
      biotype = target$biotype %||% NA_character_,
      tractability = tract_txt,
      known_drug_count = known_drug_count,
      top_known_drugs = top_drugs,
      open_targets_status = "ok"
    )
  })
  saveRDS(out, cache_path)
  out
}

module_membership_scores <- function(wgcna_res) {
  if (is.null(wgcna_res$net) || !nrow(wgcna_res$modules)) {
    return(tibble::tibble(gene = character(), module = character(), module_score = numeric()))
  }
  # Intramodular membership proxy: inverse module size; detailed kME can be added if the full multiExpr object is retained.
  wgcna_res$modules |>
    dplyr::count(.data$module, name = "module_size") |>
    dplyr::right_join(wgcna_res$modules, by = "module") |>
    dplyr::mutate(module_score = 1 / sqrt(.data$module_size))
}

make_integrative_hub_table <- function(meta_res, adjusted_meta, core, microenv, transition, wgcna_res, module_meta, cfg) {
  base <- meta_res$meta |>
    dplyr::mutate(
      gene = .data$feature,
      abs_meta_z = abs(.data$z),
      heterogeneity_penalty = pmax(.data$I2, 0) / 100,
      core_flag = .data$feature %in% core$feature,
      microenvironment_flag = .data$feature %in% microenv$feature,
      transition_flag = .data$feature %in% transition$feature
    )
  adj <- adjusted_meta$meta |> dplyr::select(feature, adjusted_meta_logFC = .data$meta_logFC, adjusted_fdr = .data$fdr)
  mod <- module_membership_scores(wgcna_res) |> dplyr::select(gene, module, module_score)

  hubs <- base |>
    dplyr::left_join(adj, by = "feature") |>
    dplyr::left_join(mod, by = "gene") |>
    dplyr::mutate(
      retained_after_tme_adjustment = !is.na(.data$adjusted_fdr) & .data$adjusted_fdr < 0.05 & sign(.data$adjusted_meta_logFC) == sign(.data$meta_logFC),
      integrated_hub_score =
        scales::rescale(.data$abs_meta_z, to = c(0, 1), from = range(.data$abs_meta_z, na.rm = TRUE)) +
        0.75 * dplyr::coalesce(.data$direction_concordance, 0) +
        0.75 * as.numeric(.data$retained_after_tme_adjustment) +
        0.50 * as.numeric(.data$core_flag) +
        0.40 * as.numeric(.data$transition_flag) +
        0.25 * as.numeric(.data$microenvironment_flag) +
        0.25 * dplyr::coalesce(.data$module_score, 0) -
        0.75 * dplyr::coalesce(.data$heterogeneity_penalty, 1)
    ) |>
    dplyr::arrange(dplyr::desc(.data$integrated_hub_score))

  ann <- annotate_opentargets(head(hubs$gene, 250), cfg)
  hubs <- hubs |> dplyr::left_join(ann, by = c("gene" = "symbol"))
  write_tsv(hubs, file.path(cfg$project$output_dir, "tables", "integrative_hub_gene_prioritization.tsv"))
  hubs
}
