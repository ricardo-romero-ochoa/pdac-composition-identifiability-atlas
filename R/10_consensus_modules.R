# Consensus modules / WGCNA ---------------------------------------------------

prepare_wgcna_multiExpr <- function(dataset_list, cfg) {
  prepared_list <- purrr::imap(dataset_list, function(ds, dsid) {
    prepared <- prepare_dataset_for_downstream_scores(ds)
    ds$expr <- prepared$expr
    ds$meta <- prepared$meta
    ds
  })
  eligible <- prepared_list[vapply(prepared_list, function(ds) ncol(ds$expr) >= (cfg$analysis$min_samples_for_wgcna %||% 40), logical(1))]
  if (length(eligible) < 2) {
    warning("Fewer than two eligible datasets for consensus WGCNA.")
    return(NULL)
  }
  genes <- Reduce(intersect, lapply(eligible, function(ds) rownames(ds$expr)))
  if (length(genes) > (cfg$analysis$top_variable_genes_for_wgcna %||% 6000)) {
    vars <- purrr::map(eligible, function(ds) matrixStats::rowVars(ds$expr[genes, , drop = FALSE], na.rm = TRUE))
    var_rank <- rowMeans(do.call(cbind, vars), na.rm = TRUE)
    names(var_rank) <- genes
    genes <- names(sort(var_rank, decreasing = TRUE))[seq_len(cfg$analysis$top_variable_genes_for_wgcna %||% 6000)]
  }
  multiExpr <- purrr::map(eligible, function(ds) {
    dat <- t(ds$expr[genes, , drop = FALSE])
    list(data = as.data.frame(dat))
  })
  list(multiExpr = multiExpr, eligible = eligible, genes = genes)
}

run_consensus_wgcna <- function(dataset_list, cfg) {
  if (!quiet_require("WGCNA")) {
    warning("WGCNA is not installed; module analysis skipped.")
    return(list(modules = tibble::tibble(), eigengenes = list(), module_de = tibble::tibble()))
  }
  WGCNA::allowWGCNAThreads()
  prep <- prepare_wgcna_multiExpr(dataset_list, cfg)
  if (is.null(prep)) return(list(modules = tibble::tibble(), eigengenes = list(), module_de = tibble::tibble()))

  gsg <- WGCNA::goodSamplesGenesMS(prep$multiExpr, verbose = 2)
  if (!gsg$allOK) {
    prep$multiExpr <- purrr::map(prep$multiExpr, function(x) list(data = x$data[gsg$goodSamples[[1]], gsg$goodGenes, drop = FALSE]))
    prep$genes <- prep$genes[gsg$goodGenes]
  }

  net <- WGCNA::blockwiseConsensusModules(
    prep$multiExpr,
    power = cfg$analysis$wgcna_power %||% 6,
    minModuleSize = cfg$analysis$wgcna_min_module_size %||% 30,
    mergeCutHeight = cfg$analysis$wgcna_merge_cut_height %||% 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 2
  )

  module_colors <- WGCNA::labels2colors(net$colors)
  modules <- tibble::tibble(gene = prep$genes, module = module_colors)
  write_tsv(modules, file.path(cfg$project$output_dir, "tables", "consensus_wgcna_gene_modules.tsv"))

  eigengenes <- purrr::imap(prep$eligible, function(ds, dsid) {
    genes <- intersect(modules$gene, rownames(ds$expr))
    dat <- t(ds$expr[genes, , drop = FALSE])
    colors <- modules$module[match(colnames(dat), modules$gene)]
    me <- WGCNA::moduleEigengenes(dat, colors = colors)$eigengenes
    me <- as.data.frame(me)
    rownames(me) <- rownames(dat)
    t(as.matrix(me))
  })

  module_de <- purrr::imap(eigengenes, function(me, dsid) {
    ds <- prep$eligible[[dsid]]
    meta <- ds$meta |> dplyr::filter(.data$sample %in% colnames(me))
    me <- me[, meta$sample, drop = FALSE]
    run_limma_contrast(me, meta, contrast = "tumor_vs_control", paired_design = ds$config$paired_design %||% FALSE) |>
      dplyr::rename(module = gene) |>
      dplyr::mutate(dataset = dsid, .before = 1)
  })
  all_module_de <- dplyr::bind_rows(module_de)
  write_tsv(all_module_de, file.path(cfg$project$output_dir, "tables", "consensus_module_de_by_dataset.tsv"))

  list(modules = modules, eigengenes = eigengenes, module_de = module_de, net = net)
}

annotate_modules <- function(wgcna_res, core, microenv, transition, cfg) {
  modules <- wgcna_res$modules
  if (!nrow(modules)) return(tibble::tibble())
  prog_lists <- list(
    core_tumor = core$feature %||% character(0),
    microenvironment = microenv$feature %||% character(0),
    transition = transition$feature %||% character(0)
  )
  universe <- unique(modules$gene)
  out <- purrr::map_dfr(unique(modules$module), function(m) {
    genes <- modules$gene[modules$module == m]
    purrr::map_dfr(names(prog_lists), function(pn) {
      hits <- intersect(genes, prog_lists[[pn]])
      bg_hits <- intersect(universe, prog_lists[[pn]])
      p <- stats::phyper(length(hits) - 1, length(bg_hits), length(universe) - length(bg_hits), length(genes), lower.tail = FALSE)
      tibble::tibble(module = m, annotation = pn, n_module = length(genes), n_hits = length(hits), p = p, genes = paste(hits, collapse = ";"))
    })
  }) |>
    dplyr::mutate(fdr = bh(.data$p)) |>
    dplyr::arrange(.data$fdr)
  write_tsv(out, file.path(cfg$project$output_dir, "tables", "consensus_module_program_enrichment.tsv"))
  out
}

meta_analyze_modules <- function(wgcna_res, cfg) {
  if (!length(wgcna_res$module_de)) return(NULL)
  meta_analyze_de(wgcna_res$module_de, cfg, feature_col = "module", out_prefix = "module")
}
