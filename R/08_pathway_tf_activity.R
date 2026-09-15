# Pathway and regulator activity ---------------------------------------------

get_hallmark_sets <- function() {
  # Do not use split(.$gs_name) after the base pipe: the magrittr dot pronoun is
  # magrittr-specific and fails inside targets workers with "object '.' not found".
  msig <- tryCatch({
    msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }, warning = function(w) {
    suppressWarnings(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  })
  required <- c("gs_name", "gene_symbol")
  missing <- setdiff(required, colnames(msig))
  if (length(missing)) {
    stop("msigdbr Hallmark output is missing column(s): ", paste(missing, collapse = ", "),
         ". Check the installed msigdbr version.", call. = FALSE)
  }
  msig <- msig |>
    dplyr::select("gs_name", "gene_symbol") |>
    dplyr::mutate(gene_symbol = toupper(.data$gene_symbol)) |>
    dplyr::filter(!is.na(.data$gs_name), !is.na(.data$gene_symbol), nzchar(.data$gene_symbol))
  split(msig$gene_symbol, msig$gs_name) |>
    purrr::map(unique)
}

run_gsva_one <- function(ds) {
  sets <- get_hallmark_sets()
  sets <- purrr::map(sets, ~ intersect(.x, rownames(ds$expr)))
  sets <- sets[lengths(sets) >= 10]
  # Use modern GSVA interface when available; fall back to legacy call.
  res <- tryCatch({
    par <- GSVA::ssgseaParam(ds$expr, sets, normalize = TRUE)
    GSVA::gsva(par)
  }, error = function(e) {
    GSVA::gsva(ds$expr, sets, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
  })
  as.matrix(res)
}

normalise_activity_result <- function(res, sample_ids, dataset = NA_character_, collection = NA_character_) {
  # Convert a method-specific activity object to the common convention used by
  # the rest of the pipeline: features/pathways/regulators x samples.
  #
  # This is intentionally conservative. Some PROGENy versions interpret input
  # matrices differently and may return pathways x genes instead of pathways x
  # samples when expression is transposed incorrectly. Such objects are rejected
  # here instead of being sent to limma with empty metadata.
  if (is.null(res)) return(NULL)

  if (is.data.frame(res)) {
    char_cols <- names(res)[vapply(res, is.character, logical(1))]
    num_cols <- names(res)[vapply(res, is.numeric, logical(1))]
    if (length(char_cols) >= 1L && length(num_cols) >= 1L && is.null(rownames(res))) {
      id_col <- char_cols[[1]]
      ids <- as.character(res[[id_col]])
      if (!anyDuplicated(ids) && all(nzchar(ids))) {
        tmp <- as.data.frame(res[, num_cols, drop = FALSE])
        rownames(tmp) <- ids
        res <- tmp
      }
    }
  }

  mat <- suppressWarnings(as.matrix(res))
  storage.mode(mat) <- "numeric"

  if (is.null(rownames(mat)) || is.null(colnames(mat))) return(NULL)

  sample_ids <- as.character(sample_ids)
  col_overlap <- length(intersect(colnames(mat), sample_ids))
  row_overlap <- length(intersect(rownames(mat), sample_ids))

  if (col_overlap == 0L && row_overlap == 0L) return(NULL)

  # If samples are rows, transpose. If both overlap, choose the larger overlap.
  if (row_overlap > col_overlap) mat <- t(mat)

  common <- intersect(colnames(mat), sample_ids)
  if (length(common) < 3L) return(NULL)

  mat <- mat[, common, drop = FALSE]
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 1L) return(NULL)

  mat
}

run_progeny_one <- function(ds) {
  if (!quiet_require("progeny")) return(NULL)

  sample_ids <- colnames(ds$expr)

  # PROGENy input orientation has differed across wrappers/examples. In this
  # repo, ds$expr is genes x samples. Try that first because it returns
  # samples x pathways in the commonly installed progeny API. If it cannot be
  # aligned to sample IDs, try the transposed matrix.
  attempts <- list(genes_by_samples = ds$expr, samples_by_genes = t(ds$expr))

  results <- purrr::imap(attempts, function(expr_in, nm) {
    res <- tryCatch({
      progeny::progeny(expr_in, scale = TRUE, organism = "Human", top = 500, perm = 1)
    }, error = function(e) {
      warning(ds$dataset, ": PROGENy failed with ", nm, " input: ", conditionMessage(e))
      NULL
    })
    mat <- normalise_activity_result(res, sample_ids, dataset = ds$dataset, collection = paste0("progeny:", nm))
    if (is.null(mat)) return(NULL)
    attr(mat, "progeny_input_orientation") <- nm
    mat
  })

  results <- purrr::compact(results)
  if (!length(results)) {
    warning(
      ds$dataset,
      ": PROGENy returned no matrix that could be aligned to sample IDs. ",
      "The collection will be skipped for this dataset."
    )
    return(NULL)
  }

  # Prefer the result with the largest number of matched samples; ties usually
  # favour the first attempt, genes_by_samples.
  n_samples <- vapply(results, ncol, integer(1))
  results[[which.max(n_samples)]]
}

run_dorothea_one <- function(ds) {
  if (!quiet_require("dorothea") || !quiet_require("viper")) return(NULL)
  data("dorothea_hs", package = "dorothea", envir = environment())
  regulon_df <- dorothea_hs |>
    dplyr::filter(.data$confidence %in% c("A", "B", "C")) |>
    dplyr::mutate(target = toupper(.data$target), tf = toupper(.data$tf))
  regulon <- dorothea::df2regulon(regulon_df)
  res <- tryCatch({
    viper::viper(ds$expr, regulon, verbose = FALSE)
  }, error = function(e) {
    warning(ds$dataset, ": DoRothEA/VIPER failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(res)) return(NULL)
  as.matrix(res)
}

run_activity_one <- function(ds) {
  prepared <- prepare_dataset_for_downstream_scores(ds)
  ds_prepared <- ds
  ds_prepared$expr <- prepared$expr
  ds_prepared$meta <- prepared$meta
  list(
    hallmark_ssgsea = run_gsva_one(ds_prepared),
    progeny = run_progeny_one(ds_prepared),
    dorothea_viper = run_dorothea_one(ds_prepared)
  )
}

run_all_activity <- function(dataset_list, cfg) {
  out <- purrr::map(dataset_list, run_activity_one)
  saveRDS(out, file.path(cfg$project$output_dir, "objects", "activity_scores_by_dataset.rds"))
  out
}

orient_activity_matrix <- function(mat, meta, dataset = NA_character_, collection = NA_character_) {
  # Activity engines do not all guarantee the same orientation across versions.
  # GSVA usually returns features x samples, but PROGENy/VIPER wrappers and
  # fallback calls may return samples x features. limma needs features x samples.
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(
      dataset, ": activity matrix for ", collection,
      " lacks rownames or colnames; cannot align it to metadata samples.",
      call. = FALSE
    )
  }
  sample_ids <- as.character(meta$sample)
  n_col_overlap <- length(intersect(colnames(mat), sample_ids))
  n_row_overlap <- length(intersect(rownames(mat), sample_ids))

  if (n_col_overlap == 0L && n_row_overlap == 0L) {
    stop(
      dataset, ": activity matrix for ", collection,
      " has no overlap with metadata sample IDs. ",
      "First activity columns: ", paste(utils::head(colnames(mat), 5), collapse = ", "),
      "; first activity rows: ", paste(utils::head(rownames(mat), 5), collapse = ", "),
      "; first metadata samples: ", paste(utils::head(sample_ids, 5), collapse = ", "),
      call. = FALSE
    )
  }

  if (n_row_overlap > n_col_overlap) {
    mat <- t(mat)
  }

  common <- intersect(colnames(mat), sample_ids)
  if (length(common) < 3L) {
    stop(
      dataset, ": fewer than three samples overlap between activity matrix and metadata for ",
      collection, " after orientation check; overlap=", length(common), ".",
      call. = FALSE
    )
  }

  meta2 <- meta |>
    dplyr::filter(.data$sample %in% common) |>
    dplyr::arrange(match(.data$sample, common))
  mat2 <- mat[, meta2$sample, drop = FALSE]

  if (anyDuplicated(meta2$sample) || anyDuplicated(colnames(mat2))) {
    stop(dataset, ": duplicate sample IDs remain after activity alignment for ", collection, call. = FALSE)
  }
  list(mat = mat2, meta = meta2)
}

run_activity_de <- function(activity_list, dataset_list, cfg, collection = "hallmark_ssgsea", contrast = "tumor_vs_control") {
  out <- purrr::imap(activity_list, function(a, dsid) {
    mat <- a[[collection]]
    if (is.null(mat)) return(NULL)
    ds <- dataset_list[[dsid]]
    prepared <- prepare_dataset_for_downstream_scores(ds)
    aligned <- orient_activity_matrix(mat, prepared$meta, dataset = dsid, collection = collection)
    run_limma_contrast(
      aligned$mat,
      aligned$meta,
      contrast = contrast,
      paired_design = ds$config$paired_design %||% FALSE
    ) |>
      dplyr::rename(term = gene) |>
      dplyr::mutate(dataset = dsid, collection = collection, .before = 1)
  })
  out <- purrr::compact(out)
  all <- dplyr::bind_rows(out)
  write_tsv(all, file.path(cfg$project$output_dir, "tables", paste0("activity_de_", collection, ".tsv")))
  out
}

run_all_activity_de_and_meta <- function(activity_list, dataset_list, cfg) {
  collections <- c("hallmark_ssgsea", "progeny", "dorothea_viper")
  out <- purrr::map(collections, function(coll) {
    de <- run_activity_de(activity_list, dataset_list, cfg, collection = coll)
    if (!length(de)) return(NULL)
    names(de) <- unique(dplyr::bind_rows(de)$dataset)
    meta <- meta_analyze_de(de, cfg, feature_col = "term", out_prefix = paste0("activity_", coll))
    list(de = de, meta = meta)
  })
  names(out) <- collections
  purrr::compact(out)
}
