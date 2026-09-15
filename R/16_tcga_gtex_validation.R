# TCGA-PAAD / GTEx pancreas validation layer.
# Preferred input: preprocessed local expression matrix (genes x samples) and
# metadata table. Optional UCSC Xena Toil downloader is included for convenience.
#
# v3.5 note: the UCSC Xena Toil expression matrix is large. This module now
# performs strict PAAD/pancreas sample selection and streams only the genes needed
# for the atlas signatures, preventing OS-level "Terminated" kills on laptops.

get_tcga_gtex_cfg <- function(cfg) get_validation_subcfg(cfg, "tcga_gtex")

.download_if_needed <- function(url, dest, force = FALSE) {
  if (!file.exists(dest) || force) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    message("Downloading ", url, " -> ", dest)
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  }
  dest
}

map_ensembl_ids_to_symbols_safe <- function(ids) {
  ids <- safe_to_utf8(as.character(ids))
  ens <- sub("\\..*$", "", ids)
  out <- ids
  looks_ensembl <- grepl("^ENSG[0-9]+$", ens)
  if (!any(looks_ensembl)) return(out)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(out)
  }
  keys <- unique(ens[looks_ensembl])
  mp <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = keys,
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    )),
    error = function(e) NULL
  )
  if (is.null(mp) || !nrow(mp) || !all(c("ENSEMBL", "SYMBOL") %in% names(mp))) return(out)
  mp <- mp[!is.na(mp$SYMBOL) & nzchar(mp$SYMBOL), c("ENSEMBL", "SYMBOL"), drop = FALSE]
  mp <- mp[!duplicated(mp$ENSEMBL), , drop = FALSE]
  sym <- mp$SYMBOL[match(ens, mp$ENSEMBL)]
  replace <- !is.na(sym) & nzchar(sym)
  out[replace] <- sym[replace]
  out
}

normalise_xena_gene_ids <- function(x) {
  x <- safe_to_utf8(as.character(x))
  maybe_mapped <- map_ensembl_ids_to_symbols_safe(x)
  normalise_gene_symbols(maybe_mapped)
}

signature_gene_universe <- function(catalog) {
  modules <- names(catalog)[!names(catalog) %in% "tabs"]
  genes <- unlist(lapply(modules, function(nm) {
    sig <- catalog[[nm]]
    c(sig$up %||% character(), sig$down %||% character(), names(sig$weights %||% numeric()))
  }), use.names = FALSE)
  genes <- normalise_gene_symbols(genes)
  unique(genes[!is.na(genes) & nzchar(genes)])
}

symbol_universe_to_ensembl <- function(symbols) {
  symbols <- normalise_gene_symbols(symbols)
  if (!length(symbols)) return(character())
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(character())
  }
  mp <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(symbols),
      keytype = "SYMBOL",
      columns = c("ENSEMBL")
    )),
    error = function(e) NULL
  )
  if (is.null(mp) || !nrow(mp) || !"ENSEMBL" %in% names(mp)) return(character())
  unique(mp$ENSEMBL[!is.na(mp$ENSEMBL) & nzchar(mp$ENSEMBL)])
}

standardize_xena_gene_matrix <- function(expr, gene_col = NULL) {
  expr <- as.data.frame(expr, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(gene_col)) {
    gene_col <- intersect(c("sample", "gene", "Gene", "Name", "id"), names(expr))[1]
    if (is.na(gene_col)) gene_col <- names(expr)[1]
  }
  expr[[gene_col]] <- normalise_xena_gene_ids(expr[[gene_col]])
  standardize_gene_matrix(expr, gene_col = gene_col)
}

infer_tcga_gtex_groups <- function(pheno) {
  pheno <- clean_character_columns(pheno)
  names(pheno) <- make.names(names(pheno), unique = TRUE)
  sample_col <- intersect(c("sample", "Sample", "sampleID", "sample_id", "X_sample", "X.sample", "id"), names(pheno))[1]
  if (is.na(sample_col)) sample_col <- names(pheno)[1]
  pheno$sample <- safe_to_utf8(pheno[[sample_col]])
  collapsed <- collapse_rows_lower(pheno)

  is_tcga <- grepl("tcga", collapsed, ignore.case = TRUE) | grepl("^TCGA", pheno$sample, ignore.case = TRUE)
  is_gtex <- grepl("gtex", collapsed, ignore.case = TRUE) | grepl("^GTEX", pheno$sample, ignore.case = TRUE)
  is_pancreas <- grepl("pancre", collapsed)
  is_paad <- grepl("paad|pancreatic adenocarcinoma|pancreatic ductal adenocarcinoma|pancreas", collapsed)

  # TCGA primary tumor barcodes contain -01 after the patient barcode. In Xena,
  # sample type is also usually represented in phenotype fields.
  is_primary_tumor <- grepl("primary tumor|primary solid tumor|tumou?r|pdac|paad", collapsed) |
    grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", pheno$sample)

  pheno$condition <- NA_character_
  pheno$cohort <- NA_character_

  # Strictly restrict TCGA to PAAD/pancreas rows. Earlier versions had a broad
  # fallback that could select every TCGA primary tumor across cancer types, which
  # made the Xena expression read enormous and could trigger an OS "Terminated".
  tcga_paad_tumor <- is_tcga & is_paad & is_primary_tumor
  gtex_pancreas <- is_gtex & is_pancreas

  pheno$condition[tcga_paad_tumor] <- "tumor"
  pheno$cohort[tcga_paad_tumor] <- "TCGA_PAAD"
  pheno$condition[gtex_pancreas] <- "control"
  pheno$cohort[gtex_pancreas] <- "GTEx_pancreas"

  out <- pheno |>
    dplyr::filter(!is.na(.data$condition)) |>
    dplyr::distinct(.data$sample, .keep_all = TRUE)

  if (!nrow(out)) {
    stop(
      "No TCGA-PAAD / GTEx pancreas samples could be inferred from the Xena phenotype file. ",
      "Inspect data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz and, if needed, ",
      "provide preprocessed metadata with columns sample, condition, cohort.",
      call. = FALSE
    )
  }
  msg <- out |>
    dplyr::count(.data$cohort, .data$condition, name = "n") |>
    dplyr::mutate(txt = paste0(.data$cohort, "/", .data$condition, "=", .data$n)) |>
    dplyr::pull(.data$txt) |>
    paste(collapse = "; ")
  message("Selected TCGA/GTEx samples: ", msg)
  out
}

read_xena_header <- function(expr_file) {
  con <- if (grepl("\\.gz$", expr_file, ignore.case = TRUE)) gzfile(expr_file, open = "rt") else file(expr_file, open = "rt")
  on.exit(close(con), add = TRUE)
  h <- readLines(con, n = 1, warn = FALSE)
  if (!length(h)) stop("Expression file appears empty: ", expr_file, call. = FALSE)
  strsplit(h, "\t", fixed = TRUE)[[1]]
}

read_xena_selected_matrix_stream <- function(expr_file, sample_ids, gene_symbols, gene_col_candidates = c("sample", "gene", "Gene", "Name", "id"), chunk_size = 1000) {
  if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file, call. = FALSE)
  header <- read_xena_header(expr_file)
  gene_col <- intersect(gene_col_candidates, header)[1]
  if (is.na(gene_col)) gene_col <- header[1]
  gene_idx <- match(gene_col, header)
  selected <- intersect(unique(as.character(sample_ids)), header)
  if (length(selected) < 10) {
    stop("Too few selected TCGA/GTEx samples were found in expression matrix header: ", length(selected), call. = FALSE)
  }
  sample_idx <- match(selected, header)

  wanted_symbols <- normalise_gene_symbols(gene_symbols)
  wanted_symbols <- unique(wanted_symbols[!is.na(wanted_symbols) & nzchar(wanted_symbols)])
  wanted_ensembl <- symbol_universe_to_ensembl(wanted_symbols)
  wanted <- unique(c(wanted_symbols, wanted_ensembl))
  if (length(wanted) < 3) {
    stop("The atlas signature gene universe is empty or could not be built; cannot safely subset the Xena matrix.", call. = FALSE)
  }

  message(
    "Streaming Xena expression subset: ", length(selected), " samples, ",
    length(wanted_symbols), " requested atlas symbols."
  )

  con <- if (grepl("\\.gz$", expr_file, ignore.case = TRUE)) gzfile(expr_file, open = "rt") else file(expr_file, open = "rt")
  on.exit(close(con), add = TRUE)
  invisible(readLines(con, n = 1, warn = FALSE))

  idx <- c(gene_idx, sample_idx)
  rows <- list()
  total <- 0L
  kept <- 0L
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    total <- total + length(lines)
    # Xena gene column is normally the first column; this fast path avoids
    # splitting every large row unless the gene is needed.
    if (identical(gene_idx, 1L)) {
      raw_gene <- sub("\t.*$", "", lines)
    } else {
      first_split <- strsplit(lines, "\t", fixed = TRUE)
      raw_gene <- vapply(first_split, function(z) z[[gene_idx]], character(1))
    }
    gene_symbol_key <- normalise_gene_symbols(raw_gene)
    gene_ens_key <- sub("\\..*$", "", safe_to_utf8(raw_gene))
    keep <- gene_symbol_key %in% wanted | gene_ens_key %in% wanted
    if (any(keep)) {
      pieces <- strsplit(lines[keep], "\t", fixed = TRUE)
      add <- lapply(pieces, function(z) {
        if (length(z) < max(idx)) return(NULL)
        z[idx]
      })
      add <- add[!vapply(add, is.null, logical(1))]
      if (length(add)) {
        rows <- c(rows, add)
        kept <- kept + length(add)
      }
    }
  }

  if (!length(rows)) {
    stop(
      "No atlas signature genes were found in the Xena expression matrix. ",
      "The first column may use identifiers that could not be mapped. ",
      "Install org.Hs.eg.db or provide a preprocessed expression file with gene symbols.",
      call. = FALSE
    )
  }

  mat_chr <- do.call(rbind, rows)
  out <- as.data.frame(mat_chr, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(out) <- c(gene_col, selected)
  message("Read ", total, " Xena rows and retained ", kept, " atlas-relevant rows before duplicate-gene collapse.")
  standardize_xena_gene_matrix(out, gene_col = gene_col)
}

read_xena_selected_matrix <- function(expr_file, sample_ids, gene_symbols = NULL, gene_col_candidates = c("sample", "gene", "Gene", "Name", "id"), memory_safe = TRUE) {
  if (isTRUE(memory_safe) && !is.null(gene_symbols)) {
    return(read_xena_selected_matrix_stream(expr_file, sample_ids, gene_symbols = gene_symbols, gene_col_candidates = gene_col_candidates))
  }
  if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file, call. = FALSE)
  header <- names(data.table::fread(expr_file, nrows = 0, data.table = FALSE))
  gene_col <- intersect(gene_col_candidates, header)[1]
  if (is.na(gene_col)) gene_col <- header[1]
  selected <- intersect(sample_ids, header)
  if (length(selected) < 10) {
    stop("Too few selected TCGA/GTEx samples were found in expression matrix header: ", length(selected), call. = FALSE)
  }
  data.table::fread(expr_file, select = c(gene_col, selected), data.table = FALSE) |>
    standardize_xena_gene_matrix(gene_col = gene_col)
}

load_tcga_gtex_preprocessed <- function(cfg) {
  vcfg <- get_tcga_gtex_cfg(cfg)
  expr_file <- vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"
  meta_file <- vcfg$metadata_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv"
  expr <- read_table_auto(expr_file) |> standardize_gene_matrix()
  meta <- read_table_auto(meta_file, required = c("sample", "condition")) |>
    dplyr::filter(.data$sample %in% colnames(expr)) |>
    dplyr::arrange(match(.data$sample, colnames(expr)))
  expr <- expr[, meta$sample, drop = FALSE]
  list(expr = expr, meta = meta)
}

prepare_tcga_gtex_from_xena <- function(cfg) {
  vcfg <- get_tcga_gtex_cfg(cfg)
  force <- isTRUE(vcfg$force_download %||% FALSE)
  expr_url <- vcfg$xena_expression_url %||% "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz"
  pheno_url <- vcfg$xena_phenotype_url %||% "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"
  expr_file <- vcfg$xena_expression_file %||% "data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz"
  pheno_file <- vcfg$xena_phenotype_file %||% "data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz"
  if (isTRUE(vcfg$download %||% FALSE)) {
    .download_if_needed(pheno_url, pheno_file, force = force)
    .download_if_needed(expr_url, expr_file, force = force)
  }
  pheno <- data.table::fread(pheno_file, data.table = FALSE)
  pheno <- clean_character_columns(pheno)
  meta <- infer_tcga_gtex_groups(pheno)
  max_samples <- as.integer(vcfg$max_xena_samples %||% 1000L)
  if (is.finite(max_samples) && nrow(meta) > max_samples) {
    stop(
      "TCGA/GTEx sample inference selected ", nrow(meta), " samples, which exceeds max_xena_samples=", max_samples, ". ",
      "This usually means the phenotype parser selected all TCGA cancers instead of TCGA-PAAD. ",
      "Inspect the selected metadata or provide preprocessed TCGA/GTEx files.",
      call. = FALSE
    )
  }
  catalog <- build_signature_catalog(cfg)
  genes_needed <- signature_gene_universe(catalog)
  memory_safe <- isTRUE(vcfg$memory_safe_xena %||% TRUE)
  expr <- read_xena_selected_matrix(expr_file, meta$sample, gene_symbols = genes_needed, memory_safe = memory_safe)
  meta <- meta |> dplyr::filter(.data$sample %in% colnames(expr)) |> dplyr::arrange(match(.data$sample, colnames(expr)))
  expr <- expr[, meta$sample, drop = FALSE]
  out_expr <- vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"
  out_meta <- vcfg$metadata_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv"
  write_tsv_safe(tibble::as_tibble(expr, rownames = "gene"), out_expr)
  write_tsv_safe(meta, out_meta)
  list(expr = expr, meta = meta)
}

load_or_prepare_tcga_gtex <- function(cfg) {
  vcfg <- get_tcga_gtex_cfg(cfg)
  expr_file <- vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz"
  meta_file <- vcfg$metadata_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv"
  if (file.exists(expr_file) && file.exists(meta_file)) return(load_tcga_gtex_preprocessed(cfg))

  xena_expr_file <- vcfg$xena_expression_file %||% "data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz"
  xena_pheno_file <- vcfg$xena_phenotype_file %||% "data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz"
  xena_requested <- isTRUE(vcfg$source %in% c("xena_toil", "xena"))
  xena_local_ready <- file.exists(xena_expr_file) && file.exists(xena_pheno_file)
  xena_download <- isTRUE(vcfg$download %||% FALSE)
  allow_in_pipeline <- isTRUE(vcfg$allow_in_pipeline_xena %||% FALSE)
  if (xena_requested && (xena_local_ready || xena_download)) {
    if (allow_in_pipeline) {
      return(prepare_tcga_gtex_from_xena(cfg))
    }
    if (isTRUE(vcfg$skip_if_missing %||% TRUE)) {
      message(
        "TCGA/GTEx preprocessed subset not found. Xena files are present, but ",
        "allow_in_pipeline_xena=false to avoid OS-level kills. Run:\n",
        "  Rscript scripts/prepare_tcga_gtex_subset.R\n",
        "then rerun TCGA/GTEx validation."
      )
      return(NULL)
    }
    stop(
      "TCGA/GTEx preprocessed subset not found. To avoid loading the full Xena matrix inside targets, run:\n",
      "  Rscript scripts/prepare_tcga_gtex_subset.R\n",
      "This creates expression_file + metadata_file, then run scripts/run_tcga_gtex_validation.R.",
      call. = FALSE
    )
  }

  if (isTRUE(vcfg$skip_if_missing %||% TRUE)) return(NULL)
  stop("TCGA/GTEx validation data not found. Provide expression_file + metadata_file, provide local Xena files, or set validation.tcga_gtex.download=true.", call. = FALSE)
}

write_tcga_gtex_skip_outputs <- function(cfg, reason) {
  msg <- tibble::tibble(status = "skipped", analysis = "tcga_gtex", reason = reason)
  write_tsv_safe(msg, file.path(external_table_dir(cfg), "tcga_gtex_module_scores.tsv"))
  write_tsv_safe(msg, file.path(external_table_dir(cfg), "tcga_gtex_module_validation.tsv"))
  write_tsv_safe(msg, file.path(external_table_dir(cfg), "tcga_survival_results.tsv"))
  saveRDS(list(status = "skipped", reason = reason), file.path(external_object_dir(cfg), "tcga_gtex_validation.rds"))
  msg
}

run_tcga_gtex_module_validation <- function(cfg) {
  dat <- load_or_prepare_tcga_gtex(cfg)
  if (is.null(dat)) {
    return(write_tcga_gtex_skip_outputs(
      cfg,
      "TCGA/GTEx validation data not found. Provide expression_file + metadata_file or set validation.tcga_gtex.download=true."
    ))
  }
  catalog <- build_signature_catalog(cfg)
  sigdiag <- signature_catalog_diagnostics(dat$expr, catalog)
  write_tsv_safe(sigdiag, file.path(external_table_dir(cfg), "tcga_gtex_signature_gene_overlap.tsv"))
  scores <- score_catalog_on_matrix(dat$expr, catalog)
  tests <- purrr::map_dfr(setdiff(colnames(scores), "sample"), function(m) {
    sc <- scores[[m]]; names(sc) <- scores$sample
    score_test_table(sc, dat$meta, group_col = "condition", positive = "tumor", dataset = "TCGA_PAAD_vs_GTEx_pancreas", module = m)
  }) |>
    dplyr::mutate(fdr = p.adjust(.data$p, method = "BH"))
  write_tsv_safe(scores, file.path(external_table_dir(cfg), "tcga_gtex_module_scores.tsv"))
  write_tsv_safe(tests, file.path(external_table_dir(cfg), "tcga_gtex_module_validation.tsv"))
  plot_score_boxplot(scores, dat$meta, out_file = file.path(external_figure_dir(cfg), "tcga_gtex_module_scores.png"), title = "TCGA-PAAD vs GTEx pancreas module validation")
  saveRDS(list(expr_dim = dim(dat$expr), meta = dat$meta, scores = scores, tests = tests), file.path(external_object_dir(cfg), "tcga_gtex_validation.rds"))
  tests
}

cox_covariate_terms <- function(df) {
  # v3.14: use clinical covariates when the survival file supplies them. The
  # minimal generated TCGA-PAAD file has age, stage and grade; optional columns
  # such as purity, Moffitt subtype, margin status or adjuvant therapy are used
  # automatically if present and sufficiently populated.
  terms <- character()
  if ("age" %in% names(df)) {
    df$age <- suppressWarnings(as.numeric(df$age))
    if (sum(is.finite(df$age)) >= 30) terms <- c(terms, "age")
  }
  categorical <- intersect(c("stage", "grade", "moffitt_subtype", "margin_status", "adjuvant_therapy"), names(df))
  for (nm in categorical) {
    vals <- as.character(df[[nm]])
    vals[!nzchar(vals) | is.na(vals) | vals %in% c("NA", "Not Available", "not reported", "[Not Available]")] <- NA_character_
    if (sum(!is.na(vals)) >= 30 && dplyr::n_distinct(vals[!is.na(vals)]) >= 2) {
      df[[nm]] <- factor(vals)
      terms <- c(terms, nm)
    }
  }
  numeric_optional <- intersect(c("purity", "estimate_purity", "ABSOLUTE_purity", "absolute_purity", "tumor_purity"), names(df))
  for (nm in numeric_optional) {
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    if (sum(is.finite(df[[nm]])) >= 30) terms <- c(terms, nm)
  }
  list(df = df, terms = unique(terms))
}

fit_module_cox <- function(df, module_name, model_name, covariates = character()) {
  df$module_score <- suppressWarnings(as.numeric(df[[module_name]]))
  keep_cols <- c("time", "event", "module_score", covariates)
  df <- df[, intersect(keep_cols, names(df)), drop = FALSE]
  df <- df |> dplyr::filter(is.finite(.data$module_score), is.finite(.data$time), !is.na(.data$event))
  if (length(covariates)) df <- df[stats::complete.cases(df[, c("module_score", covariates), drop = FALSE]), , drop = FALSE]
  if (nrow(df) < 30 || length(unique(df$event)) < 2) {
    return(tibble::tibble(module = module_name, model = model_name, n = nrow(df), events = sum(df$event == 1, na.rm = TRUE), covariates = paste(covariates, collapse = ";"), hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, ph_global_p = NA_real_, status = "skipped_insufficient_complete_cases"))
  }
  rhs <- paste(c("scale(module_score)", covariates), collapse = " + ")
  f <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", rhs))
  fit <- tryCatch(survival::coxph(f, data = df), error = function(e) e)
  if (inherits(fit, "error")) {
    return(tibble::tibble(module = module_name, model = model_name, n = nrow(df), events = sum(df$event == 1, na.rm = TRUE), covariates = paste(covariates, collapse = ";"), hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, ph_global_p = NA_real_, status = paste0("cox_failed: ", conditionMessage(fit))))
  }
  s <- summary(fit)
  ph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  ph_global_p <- NA_real_
  if (!is.null(ph) && !is.null(ph$table) && "GLOBAL" %in% rownames(ph$table)) ph_global_p <- unname(ph$table["GLOBAL", "p"])
  tibble::tibble(
    module = module_name,
    model = model_name,
    n = nrow(df),
    events = sum(df$event == 1, na.rm = TRUE),
    covariates = paste(covariates, collapse = ";"),
    hr = unname(s$coefficients[1, "exp(coef)"]),
    ci_low = unname(s$conf.int[1, "lower .95"]),
    ci_high = unname(s$conf.int[1, "upper .95"]),
    p = unname(s$coefficients[1, "Pr(>|z|)"]),
    ph_global_p = ph_global_p,
    status = "ok"
  )
}

run_tcga_survival_analysis <- function(cfg) {
  vcfg <- get_tcga_gtex_cfg(cfg)
  surv_file <- vcfg$survival_file %||% "data/external/tcga_gtex/tcga_paad_survival.tsv"
  score_file <- file.path(external_table_dir(cfg), "tcga_gtex_module_scores.tsv")
  if (!file.exists(surv_file) || !file.exists(score_file)) {
    msg <- tibble::tibble(status = "skipped", reason = "Provide data/external/tcga_gtex/tcga_paad_survival.tsv with sample, time, event, and optional covariates.")
    write_tsv_safe(msg, file.path(external_table_dir(cfg), "tcga_survival_results.tsv"))
    return(msg)
  }
  surv <- read_table_auto(surv_file, required = c("sample", "time", "event"))
  scores <- readr::read_tsv(score_file, show_col_types = FALSE)
  if (!"sample" %in% colnames(scores)) {
    msg <- tibble::tibble(status = "skipped", reason = "TCGA/GTEx module scores were skipped or unavailable; survival analysis not run.")
    write_tsv_safe(msg, file.path(external_table_dir(cfg), "tcga_survival_results.tsv"))
    return(msg)
  }
  dat <- dplyr::inner_join(surv, scores, by = "sample")
  dat$time <- suppressWarnings(as.numeric(dat$time))
  dat$event <- suppressWarnings(as.integer(dat$event))
  cv <- cox_covariate_terms(dat)
  dat <- cv$df
  clinical_covariates <- cv$terms
  modules <- setdiff(colnames(scores), "sample")
  res <- purrr::map_dfr(modules, function(m) {
    dplyr::bind_rows(
      fit_module_cox(dat, m, "univariate", character()),
      fit_module_cox(dat, m, "multivariable_available_clinical", clinical_covariates)
    )
  }) |>
    dplyr::group_by(.data$model) |>
    dplyr::mutate(fdr = p.adjust(.data$p, method = "BH")) |>
    dplyr::ungroup()
  write_tsv_safe(res, file.path(external_table_dir(cfg), "tcga_survival_results.tsv"))
  res
}

run_tcga_gtex_validation <- function(cfg) {
  module_res <- run_tcga_gtex_module_validation(cfg)
  surv_res <- run_tcga_survival_analysis(cfg)
  list(module_validation = module_res, survival = surv_res)
}
