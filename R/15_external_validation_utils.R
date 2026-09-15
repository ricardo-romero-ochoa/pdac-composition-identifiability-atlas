# External validation utility functions for the PDAC transcriptomic atlas.
# These functions deliberately accept local preprocessed files first and only
# download from public resources when explicitly requested in config. This keeps
# the base atlas reproducible and avoids accidental multi-GB downloads.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---- Validation-config helpers ---------------------------------------------
# These helpers live in the shared external-validation utility file because the
# stand-alone scripts intentionally source only the module they need. In v3.2,
# get_validation_cfg() was defined only in the TCGA/GTEx module, so
# scripts/run_scrna_mapping.R and scripts/run_depmap_opentargets.R failed when
# run independently.
get_validation_cfg <- function(cfg) {
  cfg$validation %||% list()
}

get_validation_subcfg <- function(cfg, section) {
  vcfg <- get_validation_cfg(cfg)
  vcfg[[section]] %||% list()
}


# ---- Encoding-safe text helpers ---------------------------------------------
# Public resources such as UCSC Xena phenotype files and DepMap sample metadata
# occasionally contain Latin-1 / Windows-1252 characters or malformed byte
# sequences. Plain tolower(as.character(x)) can then fail with
# "invalid multibyte string" on UTF-8 systems. These helpers normalize text before
# case-insensitive matching, while leaving numeric columns untouched.
safe_to_utf8 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  y <- tryCatch(iconv(x, from = "", to = "UTF-8", sub = "byte"), error = function(e) as.character(x))
  y[is.na(y)] <- ""
  y
}

safe_tolower <- function(x) {
  y <- safe_to_utf8(x)
  if (!is.character(y)) y <- as.character(y)
  y[is.na(y)] <- ""
  tolower(y)
}

clean_character_columns <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- make.unique(safe_to_utf8(names(df)))
  for (nm in names(df)) {
    if (is.character(df[[nm]]) || is.factor(df[[nm]])) {
      df[[nm]] <- safe_to_utf8(df[[nm]])
    }
  }
  df
}

collapse_rows_lower <- function(df, sep = " | ") {
  df <- clean_character_columns(df)
  if (!nrow(df)) return(character())
  txt <- as.data.frame(lapply(df, safe_tolower), stringsAsFactors = FALSE, check.names = FALSE)
  txt[] <- lapply(txt, function(x) { x[is.na(x)] <- ""; x })
  do.call(paste, c(txt, sep = sep))
}

# ---- Base-atlas cache bridge -------------------------------------------------
# The optional external pipeline is intentionally run with a separate targets
# store (_targets_external). It must not call targets::tar_read() inside a target.
# Instead, scripts/run_external_validation.R materializes the required base-atlas
# objects into ordinary RDS files before tar_make() starts, and the external
# targets read those files with readRDS().
external_cache_dir <- function(cfg = NULL) {
  d <- file.path("data", "external", "cache")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

base_target_cache_path <- function(target = "dataset_list", cfg = NULL) {
  file.path(external_cache_dir(cfg), paste0("base_", target, ".rds"))
}

read_base_target_from_store <- function(target = "dataset_list", store = "_targets") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("The targets package is required to read the completed base atlas target store.", call. = FALSE)
  }
  if (!dir.exists(store)) {
    stop(
      "Base targets store '", store, "' was not found. Run the base atlas first:\n",
      "  Rscript scripts/run_pipeline.R",
      call. = FALSE
    )
  }
  if ("tar_read_raw" %in% getNamespaceExports("targets")) {
    return(targets::tar_read_raw(name = target, store = store))
  }
  eval(substitute(
    targets::tar_read(TARGET, store = STORE),
    list(TARGET = as.name(target), STORE = store)
  ))
}

materialize_base_target_cache <- function(target = "dataset_list", cfg = NULL, store = "_targets", path = NULL, overwrite = TRUE) {
  path <- path %||% base_target_cache_path(target, cfg)
  if (file.exists(path) && !isTRUE(overwrite)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  obj <- read_base_target_from_store(target = target, store = store)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

require_base_target_cache_file <- function(target = "dataset_list", cfg = NULL, path = NULL) {
  path <- path %||% base_target_cache_path(target, cfg)
  if (!file.exists(path)) {
    stop(
      "External-validation cache file not found: ", path, "\n",
      "Run:\n",
      "  Rscript scripts/run_external_validation.R\n",
      "or first materialize the cache with:\n",
      "  Rscript -e \"source('R/00_utils.R'); source('R/15_external_validation_utils.R'); cfg <- read_config('config/atlas_config.yml'); materialize_base_target_cache('dataset_list', cfg)\"",
      call. = FALSE
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_base_target_cache <- function(target = "dataset_list", cfg = NULL, path = NULL) {
  readRDS(require_base_target_cache_file(target = target, cfg = cfg, path = path))
}

external_dir <- function(cfg, ...) {
  file.path(cfg$project$output_dir %||% "results", "external", ...)
}

external_table_dir <- function(cfg) {
  d <- file.path(cfg$project$output_dir %||% "results", "tables", "external")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

external_figure_dir <- function(cfg) {
  d <- file.path(cfg$project$output_dir %||% "results", "figures", "external")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

external_object_dir <- function(cfg) {
  d <- file.path(cfg$project$output_dir %||% "results", "objects", "external")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

read_table_auto <- function(path, required = NULL) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  ext <- tolower(path)
  if (grepl("\\.rds$", ext)) return(readRDS(path))
  if (grepl("\\.csv(\\.gz)?$", ext)) {
    x <- data.table::fread(path, data.table = FALSE)
  } else {
    x <- data.table::fread(path, data.table = FALSE, sep = "\t")
  }
  x <- clean_character_columns(x)
  if (!is.null(required)) {
    missing <- setdiff(required, colnames(x))
    if (length(missing)) stop("Missing column(s) in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x
}

write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(as.data.frame(x), path)
  invisible(path)
}

normalise_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- sub("^.*\\|", "", x)       # ENSG|SYMBOL -> SYMBOL
  x <- sub("\\..*$", "", x)       # SYMBOL.1 or ENSG.1 fallback
  toupper(trimws(x))
}

standardize_gene_matrix <- function(expr, gene_col = NULL, collapse_fun = c("mean", "max_median")) {
  collapse_fun <- match.arg(collapse_fun)
  if (is.data.frame(expr) || data.table::is.data.table(expr)) {
    expr <- as.data.frame(expr)
    if (is.null(gene_col)) {
      gene_col <- intersect(c("gene", "symbol", "gene_symbol", "Gene", "SYMBOL", "feature", "Name"), colnames(expr))[1]
      if (is.na(gene_col)) gene_col <- colnames(expr)[1]
    }
    genes <- normalise_gene_symbols(expr[[gene_col]])
    mat <- as.matrix(expr[, setdiff(colnames(expr), gene_col), drop = FALSE])
    suppressWarnings(storage.mode(mat) <- "numeric")
    rownames(mat) <- genes
  } else {
    mat <- as.matrix(expr)
    rownames(mat) <- normalise_gene_symbols(rownames(mat))
    suppressWarnings(storage.mode(mat) <- "numeric")
  }
  keep <- !is.na(rownames(mat)) & nzchar(rownames(mat)) & rowSums(is.finite(mat), na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (anyDuplicated(rownames(mat))) {
    split_idx <- split(seq_len(nrow(mat)), rownames(mat))
    mat <- do.call(rbind, lapply(names(split_idx), function(g) {
      sub <- mat[split_idx[[g]], , drop = FALSE]
      if (collapse_fun == "mean") {
        val <- colMeans(sub, na.rm = TRUE)
      } else {
        med <- matrixStats::rowMedians(sub, na.rm = TRUE)
        val <- sub[which.max(abs(med)), , drop = TRUE]
      }
      val
    }))
    rownames(mat) <- names(split_idx)
  }
  mat
}

zscore_rows_safe <- function(mat) {
  mat <- as.matrix(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  z <- sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
  z[!is.finite(z)] <- 0
  z
}

zscore_cols_safe <- function(mat) {
  t(zscore_rows_safe(t(mat)))
}

signed_module_score <- function(expr, up_genes = character(), down_genes = character(), weights = NULL, min_genes = 3) {
  expr <- standardize_gene_matrix(expr)
  z <- zscore_rows_safe(expr)
  up <- intersect(normalise_gene_symbols(up_genes), rownames(z))
  down <- intersect(normalise_gene_symbols(down_genes), rownames(z))
  if (!is.null(weights)) {
    weights <- weights[intersect(names(weights), rownames(z))]
    if (length(weights) >= min_genes) {
      w <- weights / sum(abs(weights), na.rm = TRUE)
      return(drop(crossprod(w, z[names(w), , drop = FALSE])))
    }
  }
  score <- rep(NA_real_, ncol(z)); names(score) <- colnames(z)
  if (length(up) >= min_genes) score <- rowMeans(t(z[up, , drop = FALSE]), na.rm = TRUE)
  if (length(down) >= min_genes) {
    down_score <- rowMeans(t(z[down, , drop = FALSE]), na.rm = TRUE)
    if (all(is.na(score))) score <- -down_score else score <- score - down_score
  }
  score
}

safe_wilcox_auc <- function(score, group, positive = "tumor") {
  df <- data.frame(score = as.numeric(score), group = as.character(group))
  df <- df[is.finite(df$score) & !is.na(df$group), , drop = FALSE]
  if (length(unique(df$group)) != 2 || !positive %in% df$group) return(NA_real_)
  pos <- df$score[df$group == positive]
  neg <- df$score[df$group != positive]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg), ties.method = "average")
  n1 <- length(pos); n0 <- length(neg)
  u <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  u / (n1 * n0)
}

score_test_table <- function(score, meta, group_col = "condition", positive = "tumor", dataset = "external", module = "module") {
  stopifnot(group_col %in% colnames(meta))
  df <- data.frame(sample = names(score), score = as.numeric(score), stringsAsFactors = FALSE)
  if (is.null(df$sample) || any(!nzchar(df$sample))) df$sample <- rownames(meta)[seq_along(score)]
  meta2 <- meta
  if (!"sample" %in% colnames(meta2)) meta2$sample <- rownames(meta2)
  df <- dplyr::left_join(df, meta2, by = "sample")
  df <- df |> dplyr::filter(!is.na(.data[[group_col]]), is.finite(.data$score))
  groups <- unique(as.character(df[[group_col]]))
  if (length(groups) < 2 || !positive %in% groups) {
    return(tibble::tibble(dataset = dataset, module = module, n_control = NA_integer_, n_tumor = NA_integer_, mean_control = NA_real_, mean_tumor = NA_real_, logFC = NA_real_, p = NA_real_, auc = NA_real_))
  }
  ref <- setdiff(groups, positive)[1]
  df2 <- df |> dplyr::filter(.data[[group_col]] %in% c(ref, positive))
  df2$.group_for_test <- df2[[group_col]]
  p <- tryCatch(stats::wilcox.test(score ~ .group_for_test, data = df2)$p.value, error = function(e) NA_real_)
  tibble::tibble(
    dataset = dataset,
    module = module,
    n_control = sum(df2[[group_col]] == ref),
    n_tumor = sum(df2[[group_col]] == positive),
    mean_control = mean(df2$score[df2[[group_col]] == ref], na.rm = TRUE),
    mean_tumor = mean(df2$score[df2[[group_col]] == positive], na.rm = TRUE),
    logFC = mean_tumor - mean_control,
    p = p,
    auc = safe_wilcox_auc(df2$score, df2[[group_col]], positive = positive)
  )
}

read_first_existing_tsv <- function(paths) {
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  for (p in paths) {
    if (file.exists(p)) {
      x <- readr::read_tsv(p, show_col_types = FALSE)
      attr(x, "source_file") <- p
      return(x)
    }
  }
  tibble::tibble()
}

load_atlas_signature_tables <- function(cfg) {
  tdir <- file.path(cfg$project$output_dir %||% "results", "tables")
  paths <- list(
    tier1 = file.path(tdir, "program_tier1_core_signature.tsv"),
    tumor_associated_retained = file.path(tdir, "program_tumor_intrinsic_retained_signature.tsv"),
    composition_covarying = c(
      file.path(tdir, "program_composition_covarying_signature.tsv"),
      file.path(tdir, "program_microenvironment_sensitive_signature.tsv")
    ),
    # v3.11: the base pipeline writes program_benign_to_pdac_transition.tsv,
    # while the external layer historically looked for transition_program.tsv.
    # Read either name, and also fall back to the full GSE91035 transition table.
    transition = c(
      file.path(tdir, "program_benign_to_pdac_transition_coding_symbol_proxy.tsv"),
      file.path(tdir, "transition_program.tsv"),
      file.path(tdir, "program_benign_to_pdac_transition.tsv"),
      file.path(tdir, "gse91035_normal_benign_pdac_transition.tsv")
    ),
    transition_validation = file.path(tdir, "transition_module_external_validation_summary.tsv"),
    strict_hubs = file.path(tdir, "integrative_hub_gene_prioritization_strict.tsv"),
    main_hubs = file.path(tdir, "main_text_hub_candidates.tsv")
  )
  out <- lapply(paths, read_first_existing_tsv)
  names(out) <- names(paths)
  out
}

signature_from_meta_table <- function(tbl, top_n = NULL) {
  if (!nrow(tbl)) return(list(up = character(), down = character(), weights = numeric()))
  feature_col <- intersect(c("feature", "gene", "symbol"), colnames(tbl))[1]
  if (is.na(feature_col)) return(list(up = character(), down = character(), weights = numeric()))
  tbl <- tbl |> dplyr::mutate(feature_sym = normalise_gene_symbols(.data[[feature_col]]))
  if (!is.null(top_n) && nrow(tbl) > top_n) {
    if ("z" %in% colnames(tbl)) tbl <- tbl |> dplyr::arrange(dplyr::desc(abs(.data$z)))
    tbl <- utils::head(tbl, top_n)
  }
  lfccol <- intersect(c("meta_logFC", "logFC", "effect"), colnames(tbl))[1]
  if (is.na(lfccol)) return(list(up = tbl$feature_sym, down = character(), weights = setNames(rep(1, nrow(tbl)), tbl$feature_sym)))
  up <- tbl$feature_sym[tbl[[lfccol]] > 0]
  down <- tbl$feature_sym[tbl[[lfccol]] < 0]
  weights <- setNames(tbl[[lfccol]], tbl$feature_sym)
  weights <- weights[is.finite(weights) & !is.na(names(weights)) & nzchar(names(weights))]
  list(up = unique(up), down = unique(down), weights = weights)
}

build_signature_catalog <- function(cfg, top_n_core = 100, top_n_tme = 100, top_n_intrinsic = 50) {
  tabs <- load_atlas_signature_tables(cfg)
  transition <- tabs$transition
  transition_up <- character(); transition_down <- character()
  if (nrow(transition)) {
    gene_col <- intersect(c("gene", "feature", "symbol"), colnames(transition))[1]
    if (!is.na(gene_col)) {
      # If the full GSE91035 transition table is used as the fallback source,
      # keep only the predefined benign-to-PDAC transition core. If the compact
      # program table is used, this filter is a no-op.
      if ("monotonic_class" %in% colnames(transition)) {
        transition <- transition |>
          dplyr::filter(.data$monotonic_class == "benign_to_pdac_transition")
      }
      transition <- transition |> dplyr::mutate(gene_sym = normalise_gene_symbols(.data[[gene_col]]))
      if ("program" %in% colnames(transition)) {
        prg <- safe_tolower(transition$program)
        transition_up <- transition$gene_sym[grepl("up", prg)]
        transition_down <- transition$gene_sym[grepl("down", prg)]
      } else {
        trend_col <- intersect(c("trend_logFC_per_stage", "logFC", "meta_logFC"), colnames(transition))[1]
        if (!is.na(trend_col)) {
          trend <- suppressWarnings(as.numeric(transition[[trend_col]]))
          transition_up <- transition$gene_sym[is.finite(trend) & trend > 0]
          transition_down <- transition$gene_sym[is.finite(trend) & trend < 0]
        }
      }
    }
  }
  transition_up <- unique(transition_up[!is.na(transition_up) & nzchar(transition_up)])
  transition_down <- unique(transition_down[!is.na(transition_down) & nzchar(transition_down)])
  list(
    tier1_core = signature_from_meta_table(tabs$tier1, top_n_core),
    tumor_associated_retained = signature_from_meta_table(tabs$tumor_associated_retained, top_n_intrinsic),
    composition_covarying = signature_from_meta_table(tabs$composition_covarying, top_n_tme),
    transition_up = list(up = transition_up, down = character(), weights = NULL),
    transition_down = list(up = transition_down, down = character(), weights = NULL),
    transition_composite = list(up = transition_up, down = transition_down, weights = NULL),
    tabs = tabs
  )
}

signature_catalog_diagnostics <- function(expr, catalog) {
  expr_genes <- rownames(standardize_gene_matrix(expr))
  modules <- names(catalog)[!names(catalog) %in% "tabs"]
  purrr::map_dfr(modules, function(nm) {
    sig <- catalog[[nm]]
    up <- unique(normalise_gene_symbols(sig$up %||% character()))
    down <- unique(normalise_gene_symbols(sig$down %||% character()))
    weights <- sig$weights %||% numeric()
    weight_genes <- unique(normalise_gene_symbols(names(weights)))
    defined <- unique(c(up, down, weight_genes))
    defined <- defined[!is.na(defined) & nzchar(defined)]
    available <- intersect(defined, expr_genes)
    src_name <- dplyr::case_when(
      nm == "tier1_core" ~ "tier1",
      nm == "tumor_associated_retained" ~ "tumor_associated_retained",
      nm == "composition_covarying" ~ "composition_covarying",
      grepl("^transition", nm) ~ "transition",
      TRUE ~ nm
    )
    src_file <- if (src_name %in% names(catalog$tabs)) attr(catalog$tabs[[src_name]], "source_file") else NA_character_
    tibble::tibble(
      module = nm,
      n_defined = length(defined),
      n_up_defined = length(up),
      n_down_defined = length(down),
      n_available = length(available),
      n_up_available = length(intersect(up, expr_genes)),
      n_down_available = length(intersect(down, expr_genes)),
      available_fraction = if (length(defined)) length(available) / length(defined) else NA_real_,
      source_file = src_file %||% NA_character_
    )
  })
}

score_catalog_on_matrix <- function(expr, catalog, min_genes = 3) {
  modules <- names(catalog)[!names(catalog) %in% "tabs"]
  out <- lapply(modules, function(nm) {
    sig <- catalog[[nm]]
    signed_module_score(expr, sig$up %||% character(), sig$down %||% character(), sig$weights %||% NULL, min_genes = min_genes)
  })
  names(out) <- modules
  as.data.frame(out, check.names = FALSE) |>
    tibble::rownames_to_column("sample") |>
    tibble::as_tibble()
}

plot_score_boxplot <- function(score_tbl, meta, group_col = "condition", out_file = NULL, title = "Module scores") {
  if (!"sample" %in% colnames(meta)) meta$sample <- rownames(meta)
  long <- score_tbl |>
    tidyr::pivot_longer(-sample, names_to = "module", values_to = "score") |>
    dplyr::left_join(meta, by = "sample") |>
    dplyr::filter(!is.na(.data[[group_col]]), is.finite(.data$score))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[group_col]], y = .data$score)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.7) +
    ggplot2::facet_wrap(~ module, scales = "free_y") +
    ggplot2::labs(x = group_col, y = "z-score module score", title = title) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = 11, height = 7, dpi = 300)
  }
  p
}
