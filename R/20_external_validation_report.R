# External validation report helpers -----------------------------------------

pandoc_minimum_available <- function(min_version = "2.8") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) return(FALSE)
  isTRUE(rmarkdown::pandoc_available(min_version))
}

html_escape_minimal <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

write_external_validation_report_fallback <- function(output = file.path("results", "pdac_external_validation_report.html"),
                                                      reason = "Pandoc >= 2.8 was not found.") {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  tdir <- file.path("results", "tables", "external")
  fdir <- file.path("results", "figures", "external")
  tables <- if (dir.exists(tdir)) list.files(tdir, pattern = "\\.tsv(\\.gz)?$", full.names = FALSE) else character(0)
  figures <- if (dir.exists(fdir)) list.files(fdir, pattern = "\\.(png|pdf|svg|jpg|jpeg)$", full.names = FALSE, ignore.case = TRUE) else character(0)
  li_tables <- if (length(tables)) paste0("<li><code>results/tables/external/", html_escape_minimal(tables), "</code></li>") else "<li>No external tables found yet.</li>"
  li_figures <- if (length(figures)) paste0("<li><code>results/figures/external/", html_escape_minimal(figures), "</code></li>") else "<li>No external figures found yet.</li>"
  html <- c(
    '<!DOCTYPE html>',
    '<html><head><meta charset="utf-8"/>',
    '<title>External validation of the cross-cohort PDAC transcriptomic atlas - report render skipped</title>',
    '<style>body{font-family:Arial,sans-serif;max-width:900px;margin:40px auto;line-height:1.5} code{background:#f5f5f5;padding:2px 4px} pre{background:#f5f5f5;padding:12px;overflow:auto} li{margin:3px 0}</style>',
    '</head><body>',
    '<h1>External validation of the cross-cohort PDAC transcriptomic atlas</h1>',
    '<h2>R Markdown rendering skipped</h2>',
    paste0('<p><strong>Reason:</strong> ', html_escape_minimal(reason), '</p>'),
    '<p>The external-validation computations can still complete without Pandoc. Use the tables and figures below for manuscript writing, or install Pandoc and rerun <code>Rscript scripts/render_external_validation_report.R</code> / <code>Rscript scripts/run_external_validation.R</code> to regenerate the full report.</p>',
    '<h2>External result tables detected</h2>',
    '<ul>', li_tables, '</ul>',
    '<h2>External result figures detected</h2>',
    '<ul>', li_figures, '</ul>',
    '<h2>Pandoc installation options</h2>',
    '<pre>conda install -c conda-forge pandoc\n# or, on Ubuntu/WSL if the repository version is new enough:\nsudo apt-get update && sudo apt-get install pandoc</pre>',
    '</body></html>'
  )
  writeLines(html, output, useBytes = TRUE)
  normalizePath(output, mustWork = TRUE)
}

render_external_validation_report <- function(input = "reports/pdac_external_validation_report.Rmd", quiet = TRUE) {
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("results", "tables", "external"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("results", "figures", "external"), recursive = TRUE, showWarnings = FALSE)
  out <- file.path("results", "pdac_external_validation_report.html")
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("rmarkdown is not installed; writing fallback external-validation report.", call. = FALSE)
    return(write_external_validation_report_fallback(out, reason = "The R package rmarkdown is not installed."))
  }
  if (!rmarkdown::pandoc_available("2.8")) {
    warning("Pandoc >= 2.8 was not found; writing fallback external-validation report instead of failing the targets pipeline.", call. = FALSE)
    return(write_external_validation_report_fallback(out, reason = "Pandoc >= 2.8 was not found."))
  }
  rmarkdown::render(
    input = input,
    output_file = "pdac_external_validation_report.html",
    output_dir = normalizePath("results", mustWork = TRUE),
    quiet = quiet,
    knit_root_dir = normalizePath(getwd(), mustWork = TRUE),
    envir = new.env(parent = globalenv())
  )
  normalizePath(out, mustWork = FALSE)
}

run_external_validation_summary <- function(cfg) {
  tdir <- external_table_dir(cfg)
  files <- list.files(tdir, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
  summary <- tibble::tibble(
    output = basename(files),
    path = files,
    exists = file.exists(files),
    size_bytes = file.info(files)$size
  )
  write_tsv_safe(summary, file.path(tdir, "external_validation_output_manifest.tsv"))
  summary
}
