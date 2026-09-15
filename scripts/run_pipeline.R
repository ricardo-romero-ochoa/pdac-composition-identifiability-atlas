#!/usr/bin/env Rscript

if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Install packages first: Rscript scripts/install_packages.R", call. = FALSE)
}

targets::tar_make()
