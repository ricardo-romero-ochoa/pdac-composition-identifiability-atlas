# Report rendering helpers ----------------------------------------------------

pandoc_minimum_available <- function(min_version = "2.8") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) return(FALSE)
  isTRUE(rmarkdown::pandoc_available(min_version))
}

write_no_pandoc_report <- function(output_dir = file.path("results"),
                                   output_file = "pdac_atlas_report.html",
                                   title = "Cross-cohort PDAC transcriptomic atlas",
                                   reason = "Pandoc >= 2.8 was not found.") {
  output_dir_abs <- normalizePath(ensure_dir(output_dir), winslash = "/", mustWork = TRUE)
  out <- file.path(output_dir_abs, output_file)
  html <- c(
    '<!DOCTYPE html>',
    '<html><head><meta charset="utf-8"/>',
    paste0('<title>', title, ' - report render skipped</title>'),
    '<style>body{font-family:Arial,sans-serif;max-width:900px;margin:40px auto;line-height:1.5} code{background:#f5f5f5;padding:2px 4px} pre{background:#f5f5f5;padding:12px;overflow:auto}</style>',
    '</head><body>',
    paste0('<h1>', title, '</h1>'),
    '<h2>Report rendering skipped</h2>',
    paste0('<p><strong>Reason:</strong> ', reason, '</p>'),
    '<p>The computational targets can still finish. Tables and figures are written under <code>results/tables</code> and <code>results/figures</code>. Install Pandoc and rerun the report-rendering script to generate the full R Markdown HTML report.</p>',
    '<pre>Rscript scripts/render_report.R</pre>',
    '<h2>Pandoc installation options</h2>',
    '<pre>conda install -c conda-forge pandoc\n# or, on Ubuntu/WSL if the repository version is new enough:\nsudo apt-get update && sudo apt-get install pandoc</pre>',
    '</body></html>'
  )
  writeLines(html, out, useBytes = TRUE)
  normalizePath(out, winslash = "/", mustWork = TRUE)
}

render_pdac_report <- function(input = "reports/pdac_atlas_report.Rmd",
                               output_dir = file.path("results"),
                               output_file = "pdac_atlas_report.html",
                               quiet = TRUE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("Package 'rmarkdown' is not installed; writing a minimal fallback report instead.", call. = FALSE)
    return(write_no_pandoc_report(output_dir, output_file, reason = "The R package rmarkdown is not installed."))
  }
  if (!rmarkdown::pandoc_available("2.8")) {
    warning("Pandoc >= 2.8 was not found; writing a minimal fallback report instead of failing the targets pipeline.", call. = FALSE)
    return(write_no_pandoc_report(output_dir, output_file, reason = "Pandoc >= 2.8 was not found."))
  }
  if (!file.exists(input)) {
    stop("Report template not found: ", input, call. = FALSE)
  }

  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  output_dir_abs <- normalizePath(ensure_dir(output_dir), winslash = "/", mustWork = TRUE)
  intermediates_dir <- ensure_dir(file.path(output_dir, "_rmarkdown_intermediates"))
  intermediates_dir_abs <- normalizePath(intermediates_dir, winslash = "/", mustWork = TRUE)

  # Do not pass file.path(output_dir, output_file) as output_file. When the input
  # is under reports/, rmarkdown interprets relative subdirectories with respect
  # to the input document directory and fails if reports/results does not exist.
  out <- rmarkdown::render(
    input = input,
    output_file = output_file,
    output_dir = output_dir_abs,
    intermediates_dir = intermediates_dir_abs,
    knit_root_dir = project_root,
    quiet = quiet,
    envir = new.env(parent = globalenv())
  )

  normalizePath(out, winslash = "/", mustWork = TRUE)
}
