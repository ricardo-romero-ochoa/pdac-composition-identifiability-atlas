#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/20_external_validation_report.R")
out <- render_external_validation_report(quiet = FALSE)
cat(out, "\n")
