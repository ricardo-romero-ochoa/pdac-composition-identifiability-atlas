# Single-cell / single-nucleus mapping of atlas signatures.
#
# v3.8 note: the loader prefers precomputed TSV tables over h5ad by default.
# This avoids reticulate/zellkonverter Python-environment ambiguity inside targets.
# Use scripts/prepare_scrna_tables.R to convert h5ad -> compact TSV tables first.

get_scrna_cfg <- function(cfg) get_validation_subcfg(cfg, "scrna")

scrna_unknown_tokens <- function() {
  c("", "unknown", "unk", "na", "n/a", "nan", "none", "null", "unspecified",
    "unassigned", "unannotated", "not_available", "not available")
}

scrna_clean_vector <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

scrna_has_useful_values <- function(x) {
  x <- scrna_clean_vector(x)
  any(!(tolower(x) %in% scrna_unknown_tokens()))
}

scrna_infer_column <- function(meta, explicit = NULL, candidates = character(), label = "column") {
  cols <- colnames(meta)
  if (!is.null(explicit) && length(explicit)) {
    explicit_vec <- unlist(explicit, use.names = FALSE)
    explicit_vec <- as.character(explicit_vec)
    explicit_vec <- explicit_vec[!is.na(explicit_vec) & nzchar(explicit_vec)]
    if (length(explicit_vec)) {
      for (wanted in explicit_vec) {
        hit <- cols[tolower(cols) == tolower(wanted)][1]
        if (!is.na(hit) && scrna_has_useful_values(meta[[hit]])) return(hit)
      }
      # A single explicitly configured column should fail loudly. A ranked list
      # can safely fall back to automatic inference, which keeps the repo usable
      # across references that do not share the same annotation schema.
      if (length(explicit_vec) == 1) {
        stop(
          "Configured scRNA ", label, " was not found or had no useful values: ", explicit_vec[1],
          "\nAvailable metadata columns:\n  ", paste(cols, collapse = ", "),
          call. = FALSE
        )
      } else {
        warning(
          "None of the configured scRNA ", label, " candidates was found with useful values: ",
          paste(explicit_vec, collapse = ", "),
          ". Falling back to automatic inference.",
          call. = FALSE
        )
      }
    }
  }
  lower <- stats::setNames(cols, tolower(cols))
  for (cand in candidates) {
    key <- tolower(cand)
    if (key %in% names(lower)) {
      hit <- unname(lower[[key]])
      if (scrna_has_useful_values(meta[[hit]])) return(hit)
    }
  }
  # Last-resort heuristic: choose an annotation-like column with a modest number
  # of distinct labels, but avoid sample/patient/barcode/QC columns.
  bad_patterns <- "barcode|cell$|sample|patient|donor|library|batch|ncount|nfeature|percent|mito|doublet|sex|age|stage|grade|score|leiden|cluster"
  candidate_cols <- cols[grepl("cell|type|annot|label|class|compartment|lineage|subtype|identity", tolower(cols)) &
                           !grepl(bad_patterns, tolower(cols))]
  for (hit in candidate_cols) {
    vals <- scrna_clean_vector(meta[[hit]])
    vals <- vals[!(tolower(vals) %in% scrna_unknown_tokens())]
    n_unique <- length(unique(vals))
    if (n_unique >= 2 && n_unique <= 100) return(hit)
  }
  NA_character_
}

scrna_cell_type_candidates <- function() {
  c(
    "cell_type", "cell_types", "celltype", "cell type", "CellType", "Cell_type",
    "cellType", "celltype_major", "cell_type_major", "major_cell_type",
    "major_celltype", "broad_cell_type", "broad_celltype", "cell_type_broad",
    "annotation", "annotations", "cell_annotation", "manual_annotation",
    "final_annotation", "author_cell_type", "author_celltype", "predicted_cell_type",
    "predicted.celltype", "predicted.id", "predicted_ID", "cell_ontology_class",
    "cell_ontology", "CellOntology", "subclass", "class", "compartment",
    "lineage", "cell.labels", "cell_label", "labels", "ident", "seurat_clusters"
  )
}

scrna_condition_candidates <- function() {
  c("condition", "condition_original", "disease", "Disease", "diagnosis", "group", "tissue",
    "sample_type", "tumor_normal", "status", "phenotype", "pathology")
}

scrna_sample_candidates <- function() {
  c("sample", "sample_id", "Sample", "SampleID", "patient", "patient_id", "donor", "donor_id",
    "orig.ident", "library", "library_id", "specimen", "case", "case_id")
}

write_scrna_metadata_diagnostics <- function(meta, cfg, prefix = "scrna_metadata") {
  out_dir <- external_table_dir(cfg)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  columns <- tibble::tibble(
    column = colnames(meta),
    n_non_missing = vapply(meta, function(x) sum(scrna_clean_vector(x) != ""), integer(1)),
    n_unique = vapply(meta, function(x) length(unique(scrna_clean_vector(x))), integer(1)),
    example_values = vapply(meta, function(x) {
      vals <- unique(scrna_clean_vector(x))
      vals <- vals[nzchar(vals)]
      paste(utils::head(vals, 8), collapse = "; ")
    }, character(1))
  )
  write_tsv_safe(columns, file.path(out_dir, paste0(prefix, "_columns.tsv")))
  invisible(columns)
}

scrna_table_files_exist <- function(expr_file, meta_file) {
  !is.null(expr_file) && !is.null(meta_file) &&
    nzchar(expr_file) && nzchar(meta_file) &&
    file.exists(expr_file) && file.exists(meta_file)
}

scrna_h5ad_error_message <- function(h5ad_file, original_error = NULL) {
  msg <- paste0(
    "Could not read h5ad input through R. The usual cause is that R/reticulate ",
    "or zellkonverter is using a different Python environment than the one where ",
    "anndata was installed.\n\n",
    "Recommended robust workflow:\n",
    "  python3 -m pip install --user anndata h5py pandas scipy numpy\n",
    "  Rscript scripts/prepare_scrna_tables.R\n",
    "  Rscript scripts/run_scrna_mapping.R\n\n",
    "The preparation script writes compact TSV files and the scRNA target will ",
    "prefer those tables over the h5ad file.\n\n",
    "Input h5ad: ", h5ad_file
  )
  if (!is.null(original_error)) msg <- paste0(msg, "\n\nOriginal error: ", conditionMessage(original_error))
  msg
}

load_scrna_reference <- function(cfg) {
  scfg <- get_scrna_cfg(cfg)
  expr_file <- scfg$expression_file %||% "data/external/scrna/scrna_expression.tsv.gz"
  meta_file <- scfg$metadata_file %||% "data/external/scrna/scrna_metadata.tsv"
  object_file <- scfg$object_file %||% NA_character_
  h5ad_file <- scfg$h5ad_file %||% NA_character_
  prefer_tables <- isTRUE(scfg$prefer_tables %||% TRUE)

  if (prefer_tables && scrna_table_files_exist(expr_file, meta_file)) {
    expr <- read_table_auto(expr_file) |> standardize_gene_matrix()
    meta <- read_table_auto(meta_file)
    if (!"cell" %in% colnames(meta)) meta$cell <- colnames(expr)
  } else if (!is.na(object_file) && file.exists(object_file)) {
    obj <- readRDS(object_file)
    if (inherits(obj, "Seurat")) {
      if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required to read Seurat references.", call. = FALSE)
      assay <- scfg$assay %||% SeuratObject::DefaultAssay(obj)
      slot <- scfg$slot %||% "data"
      expr <- SeuratObject::GetAssayData(obj, assay = assay, slot = slot)
      meta <- obj@meta.data |> tibble::rownames_to_column("cell")
      expr <- as.matrix(expr)
    } else if (inherits(obj, "SingleCellExperiment")) {
      assay <- scfg$assay %||% "logcounts"
      expr <- SummarizedExperiment::assay(obj, assay)
      meta <- as.data.frame(SummarizedExperiment::colData(obj)) |> tibble::rownames_to_column("cell")
      expr <- as.matrix(expr)
    } else if (is.list(obj) && all(c("expr", "meta") %in% names(obj))) {
      expr <- obj$expr; meta <- obj$meta
      if (!"cell" %in% colnames(meta)) meta$cell <- colnames(expr)
    } else {
      stop("Unsupported scRNA object. Provide Seurat, SingleCellExperiment, or list(expr, meta).", call. = FALSE)
    }
  } else if (!is.na(h5ad_file) && file.exists(h5ad_file)) {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
      stop(scrna_h5ad_error_message(h5ad_file), call. = FALSE)
    }
    sce <- tryCatch(zellkonverter::readH5AD(h5ad_file), error = function(e) e)
    if (inherits(sce, "error")) stop(scrna_h5ad_error_message(h5ad_file, sce), call. = FALSE)
    assay <- scfg$assay %||% if ("logcounts" %in% SummarizedExperiment::assayNames(sce)) "logcounts" else "X"
    expr <- as.matrix(SummarizedExperiment::assay(sce, assay))
    meta <- as.data.frame(SummarizedExperiment::colData(sce)) |> tibble::rownames_to_column("cell")
  } else if (scrna_table_files_exist(expr_file, meta_file)) {
    expr <- read_table_auto(expr_file) |> standardize_gene_matrix()
    meta <- read_table_auto(meta_file)
    if (!"cell" %in% colnames(meta)) meta$cell <- colnames(expr)
  } else {
    stop(
      "scRNA input files were not found. Provide either compact TSV tables:\n",
      "  ", expr_file, "\n",
      "  ", meta_file, "\n",
      "or an h5ad file at:\n",
      "  ", h5ad_file, "\n",
      "Then run:\n",
      "  Rscript scripts/prepare_scrna_tables.R\n",
      "  Rscript scripts/run_scrna_mapping.R",
      call. = FALSE
    )
  }

  expr <- standardize_gene_matrix(expr)
  meta <- meta |> dplyr::filter(.data$cell %in% colnames(expr)) |> dplyr::arrange(match(.data$cell, colnames(expr)))
  expr <- expr[, meta$cell, drop = FALSE]
  write_scrna_metadata_diagnostics(meta, cfg)

  explicit_ct <- scfg$cell_type_column %||% NULL
  cell_type_col <- scrna_infer_column(meta, explicit_ct, scrna_cell_type_candidates(), "cell_type_column")
  current_ct_bad <- (!"cell_type" %in% colnames(meta)) || !scrna_has_useful_values(meta$cell_type)
  if (!is.na(cell_type_col) && (current_ct_bad || !identical(cell_type_col, "cell_type"))) {
    meta$cell_type <- scrna_clean_vector(meta[[cell_type_col]])
  } else if (!"cell_type" %in% colnames(meta)) {
    meta$cell_type <- "unknown"
  } else {
    meta$cell_type <- scrna_clean_vector(meta$cell_type)
  }
  meta$cell_type[tolower(meta$cell_type) %in% scrna_unknown_tokens()] <- "unknown"

  explicit_condition <- scfg$condition_column %||% NULL
  condition_col <- scrna_infer_column(meta, explicit_condition, scrna_condition_candidates(), "condition_column")
  current_condition_bad <- (!"condition" %in% colnames(meta)) || !scrna_has_useful_values(meta$condition)
  if (!is.na(condition_col) && (current_condition_bad || !identical(condition_col, "condition"))) {
    meta$condition <- scrna_clean_vector(meta[[condition_col]])
  } else if (!"condition" %in% colnames(meta)) {
    meta$condition <- NA_character_
  } else {
    meta$condition <- scrna_clean_vector(meta$condition)
  }

  explicit_sample <- scfg$sample_column %||% NULL
  sample_col <- scrna_infer_column(meta, explicit_sample, scrna_sample_candidates(), "sample_column")
  if (!is.na(sample_col) && !"sample_id" %in% colnames(meta)) meta$sample_id <- scrna_clean_vector(meta[[sample_col]])

  if (isTRUE(scfg$fail_if_unknown_celltype %||% TRUE)) {
    non_unknown <- unique(meta$cell_type[meta$cell_type != "unknown"])
    if (!length(non_unknown)) {
      stop(
        "scRNA metadata did not yield usable cell-type labels; all cells are 'unknown'.\n",
        "Inspect results/tables/external/scrna_metadata_columns.tsv or run:\n",
        "  Rscript scripts/inspect_scrna_metadata.R\n",
        "Then set validation.scrna.cell_type_column in config/atlas_config.yml to the correct metadata column,\n",
        "or regenerate the compact TSV files after setting that option.",
        call. = FALSE
      )
    }
  }

  list(expr = expr, meta = meta)
}

downsample_scrna <- function(expr, meta, max_cells_per_celltype = 2000, seed = 20260531) {
  # Keep this function free of dplyr::n() inside slice_sample().
  # Recent dplyr versions require `n` in slice_sample() to be a constant,
  # so min(dplyr::n(), max_cells_per_celltype) fails even inside grouped data.
  if (is.null(max_cells_per_celltype) || !is.finite(max_cells_per_celltype)) {
    return(list(expr = expr, meta = meta))
  }
  max_cells_per_celltype <- as.integer(max_cells_per_celltype)
  if (is.na(max_cells_per_celltype) || max_cells_per_celltype <= 0) {
    return(list(expr = expr, meta = meta[0, , drop = FALSE]))
  }
  if (!"cell" %in% colnames(meta)) stop("scRNA metadata must contain a 'cell' column.", call. = FALSE)
  if (!"cell_type" %in% colnames(meta)) meta$cell_type <- "unknown"

  set.seed(seed)
  meta$cell_type <- as.character(meta$cell_type)
  meta$cell_type[is.na(meta$cell_type) | !nzchar(meta$cell_type)] <- "unknown"
  split_cells <- split(as.character(meta$cell), meta$cell_type, drop = TRUE)
  keep_cells <- unlist(lapply(split_cells, function(cells) {
    cells <- cells[!is.na(cells) & cells %in% colnames(expr)]
    if (length(cells) > max_cells_per_celltype) sample(cells, max_cells_per_celltype) else cells
  }), use.names = FALSE)
  keep_cells <- unique(keep_cells)

  meta2 <- meta[match(keep_cells, as.character(meta$cell)), , drop = FALSE]
  meta2 <- meta2[!is.na(meta2$cell), , drop = FALSE]
  list(expr = expr[, as.character(meta2$cell), drop = FALSE], meta = meta2)
}

pseudobulk_by_cell_type <- function(expr, meta, group_col = "cell_type") {
  stopifnot(group_col %in% colnames(meta))
  groups <- split(meta$cell, meta[[group_col]])
  mat <- do.call(cbind, lapply(groups, function(cells) {
    if (length(cells) == 1) expr[, cells, drop = TRUE] else rowMeans(expr[, cells, drop = FALSE], na.rm = TRUE)
  }))
  colnames(mat) <- names(groups)
  mat
}

scrna_dotplot_table <- function(expr, meta, genes, group_col = "cell_type", expression_threshold = 0) {
  genes <- intersect(normalise_gene_symbols(genes), rownames(expr))
  if (!length(genes)) return(tibble::tibble())
  groups <- split(meta$cell, meta[[group_col]])
  purrr::map_dfr(names(groups), function(g) {
    cells <- groups[[g]]
    sub <- expr[genes, cells, drop = FALSE]
    tibble::tibble(
      group = g,
      gene = genes,
      mean_expression = rowMeans(sub, na.rm = TRUE),
      fraction_expressing = rowMeans(sub > expression_threshold, na.rm = TRUE),
      n_cells = length(cells)
    )
  })
}

write_scrna_skip_outputs <- function(cfg, reason) {
  msg <- tibble::tibble(status = "skipped", analysis = "scrna_mapping", reason = reason)
  write_tsv_safe(msg, file.path(external_table_dir(cfg), "scrna_cell_module_scores.tsv.gz"))
  write_tsv_safe(msg, file.path(external_table_dir(cfg), "scrna_celltype_module_score_summary.tsv"))
  write_tsv_safe(msg, file.path(external_table_dir(cfg), "scrna_hub_gene_dotplot_table.tsv"))
  saveRDS(list(status = "skipped", reason = reason), file.path(external_object_dir(cfg), "scrna_mapping.rds"))
  msg
}

run_scrna_signature_mapping <- function(cfg) {
  scfg <- get_scrna_cfg(cfg)
  dat <- tryCatch(load_scrna_reference(cfg), error = function(e) e)
  if (inherits(dat, "error")) {
    if (isTRUE(scfg$skip_if_missing %||% TRUE)) {
      return(write_scrna_skip_outputs(cfg, conditionMessage(dat)))
    }
    stop(dat)
  }
  dat <- downsample_scrna(dat$expr, dat$meta, max_cells_per_celltype = scfg$max_cells_per_celltype %||% 2000, seed = cfg$project$random_seed %||% 20260531)
  catalog <- build_signature_catalog(cfg)
  sigdiag <- signature_catalog_diagnostics(dat$expr, catalog)
  write_tsv_safe(sigdiag, file.path(external_table_dir(cfg), "scrna_signature_gene_overlap.tsv"))
  scores <- score_catalog_on_matrix(dat$expr, catalog, min_genes = scfg$min_genes_per_signature %||% 3)
  meta <- dat$meta
  if (!"sample" %in% colnames(meta)) meta$sample <- meta$cell
  score_cell <- dplyr::left_join(scores |> dplyr::rename(cell = sample), meta, by = "cell")
  score_summary <- score_cell |>
    tidyr::pivot_longer(cols = setdiff(colnames(scores), "sample"), names_to = "module", values_to = "score") |>
    dplyr::group_by(.data$cell_type, .data$module) |>
    dplyr::summarise(n_cells = dplyr::n(), mean_score = mean(.data$score, na.rm = TRUE), median_score = stats::median(.data$score, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$module, dplyr::desc(.data$mean_score))
  hub_tbl <- load_atlas_signature_tables(cfg)$strict_hubs
  gene_col <- intersect(c("feature", "gene"), colnames(hub_tbl))[1]
  top_hubs <- if (!is.na(gene_col) && nrow(hub_tbl)) utils::head(normalise_gene_symbols(hub_tbl[[gene_col]]), scfg$top_hub_genes %||% 40) else character()
  dot <- scrna_dotplot_table(dat$expr, meta, top_hubs, group_col = "cell_type")
  write_tsv_safe(score_cell, file.path(external_table_dir(cfg), "scrna_cell_module_scores.tsv.gz"))
  write_tsv_safe(score_summary, file.path(external_table_dir(cfg), "scrna_celltype_module_score_summary.tsv"))
  write_tsv_safe(dot, file.path(external_table_dir(cfg), "scrna_hub_gene_dotplot_table.tsv"))
  p1 <- ggplot2::ggplot(score_summary, ggplot2::aes(x = cell_type, y = module, fill = mean_score)) +
    ggplot2::geom_tile() +
    ggplot2::labs(x = "Cell type", y = "Atlas module", fill = "Mean score", title = "Atlas module localization across scRNA cell types") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(external_figure_dir(cfg), "scrna_module_celltype_heatmap.png"), p1, width = 10, height = 6, dpi = 300)
  if (nrow(dot)) {
    p2 <- ggplot2::ggplot(dot, ggplot2::aes(x = group, y = gene, size = fraction_expressing, color = mean_expression)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Cell type", y = "Hub gene", size = "Fraction", color = "Mean expression", title = "Hub-gene localization across scRNA cell types") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file.path(external_figure_dir(cfg), "scrna_hub_gene_dotplot.png"), p2, width = 11, height = 8, dpi = 300)
  }
  saveRDS(list(score_summary = score_summary, dot = dot, meta_summary = dplyr::count(meta, cell_type, condition, name = "n_cells")), file.path(external_object_dir(cfg), "scrna_mapping.rds"))
  list(celltype_scores = score_summary, dotplot = dot)
}
