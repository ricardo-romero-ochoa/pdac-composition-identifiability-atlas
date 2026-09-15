# Utility functions -----------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

quiet_require <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_config <- function(path = "config/atlas_config.yml") {
  cfg <- yaml::read_yaml(path)
  ensure_dir(file.path(cfg$project$output_dir, "tables"))
  ensure_dir(file.path(cfg$project$output_dir, "figures"))
  ensure_dir(file.path(cfg$project$output_dir, "objects"))
  ensure_dir(file.path(cfg$project$data_dir, "raw"))
  ensure_dir(file.path(cfg$project$data_dir, "processed"))
  set.seed(cfg$project$random_seed %||% cfg$analysis$random_seed %||% 1)
  cfg
}

safe_file <- function(...) {
  x <- file.path(...)
  ensure_dir(dirname(x))
  x
}

save_rds <- function(x, path) {
  saveRDS(x, safe_file(path))
  invisible(path)
}

write_tsv <- function(x, path) {
  data.table::fwrite(as.data.frame(x), safe_file(path), sep = "\t")
  invisible(path)
}

log2_if_needed <- function(mat) {
  qx <- stats::quantile(as.numeric(mat), probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  # GEO processed arrays are often already log2; raw intensities are not.
  if (qx[[7]] > 100 || (qx[[6]] - qx[[1]] > 50 && qx[[2]] > 0)) {
    mat <- log2(pmax(mat, 1))
  }
  mat
}

zscore_rows <- function(mat) {
  mat <- as_numeric_matrix(mat)
  mat[!is.finite(mat)] <- NA_real_
  if (!nrow(mat) || !ncol(mat)) return(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  sd <- matrixStats::rowSds(mat, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0
  sd[!is.finite(sd) | sd == 0] <- 1
  z <- sweep(mat, 1, mu, FUN = "-")
  z <- sweep(z, 1, sd, FUN = "/")
  z[!is.finite(z)] <- 0
  z
}

sanitize_numeric_matrix <- function(mat,
                                    min_finite_per_row = 2,
                                    min_finite_per_col = 2,
                                    impute = c("row_median", "zero", "none"),
                                    drop_zero_variance_rows = TRUE,
                                    drop_zero_variance_cols = FALSE) {
  impute <- match.arg(impute)
  mat <- as_numeric_matrix(mat)
  mat[!is.finite(mat)] <- NA_real_
  if (!nrow(mat) || !ncol(mat)) return(mat[0, 0, drop = FALSE])

  keep_rows <- rowSums(is.finite(mat)) >= min_finite_per_row
  keep_cols <- colSums(is.finite(mat)) >= min_finite_per_col
  mat <- mat[keep_rows, keep_cols, drop = FALSE]
  if (!nrow(mat) || !ncol(mat)) return(mat)

  if (impute == "row_median") {
    for (i in seq_len(nrow(mat))) {
      bad <- !is.finite(mat[i, ])
      if (any(bad)) {
        med <- stats::median(mat[i, ], na.rm = TRUE)
        if (!is.finite(med)) med <- 0
        mat[i, bad] <- med
      }
    }
  } else if (impute == "zero") {
    mat[!is.finite(mat)] <- 0
  }

  if (drop_zero_variance_rows && ncol(mat) > 1) {
    rv <- matrixStats::rowVars(mat, na.rm = TRUE)
    mat <- mat[is.finite(rv) & rv > 0, , drop = FALSE]
  }
  if (drop_zero_variance_cols && nrow(mat) > 1) {
    cv <- matrixStats::colVars(mat, na.rm = TRUE)
    mat <- mat[, is.finite(cv) & cv > 0, drop = FALSE]
  }
  mat
}


normalise_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- sub("^.*\\|", "", x)       # ENSG|SYMBOL -> SYMBOL
  x <- sub("\\..*$", "", x)       # SYMBOL.1 or ENSG.1 fallback
  toupper(trimws(x))
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- gsub("///", ";", x, fixed = TRUE)
  x <- gsub("//", ";", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x[x %in% c("", "---", "NA", "na", "null", "NULL")] <- NA_character_
  # keep first listed HGNC-like symbol when probes map to multiple records
  x <- vapply(strsplit(x, ";|,| / | \\| "), function(z) trimws(z[which(nzchar(trimws(z)))[1]] %||% NA_character_), character(1))
  toupper(x)
}

collapse_by_symbol <- function(expr, symbols, strategy = c("mean", "max_iqr")) {
  strategy <- match.arg(strategy)
  symbols <- clean_symbol(symbols)
  keep <- !is.na(symbols) & nzchar(symbols)
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (!length(symbols)) stop("No valid gene symbols after probe annotation.")
  if (strategy == "mean") {
    collapsed <- limma::avereps(expr, ID = symbols)
  } else {
    split_idx <- split(seq_along(symbols), symbols)
    chosen <- vapply(split_idx, function(ii) {
      iqr <- matrixStats::rowIQRs(expr[ii, , drop = FALSE], na.rm = TRUE)
      ii[which.max(iqr)]
    }, integer(1))
    collapsed <- expr[chosen, , drop = FALSE]
    rownames(collapsed) <- names(chosen)
  }
  collapsed[order(rownames(collapsed)), , drop = FALSE]
}

first_present_col <- function(df, candidates) {
  nms <- names(df)
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[[1]] else NA_character_
}

extract_regex_group <- function(x, pattern, group = 2) {
  m <- regexec(pattern, x, ignore.case = TRUE, perl = TRUE)
  z <- regmatches(x, m)
  vapply(z, function(a) if (length(a) >= group) a[[group]] else NA_character_, character(1))
}

extract_regex_last_group <- function(x, pattern) {
  m <- regexec(pattern, x, ignore.case = TRUE, perl = TRUE)
  z <- regmatches(x, m)
  vapply(z, function(a) if (length(a) >= 2) a[[length(a)]] else NA_character_, character(1))
}

mode_or_na <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

as_numeric_matrix <- function(x) {
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

bh <- function(p) stats::p.adjust(p, method = "BH")

message_rule <- function(txt) {
  message("\n", paste(rep("=", 80), collapse = ""), "\n", txt, "\n", paste(rep("=", 80), collapse = ""))
}
