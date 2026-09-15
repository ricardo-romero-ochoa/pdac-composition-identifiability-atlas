#!/usr/bin/env Rscript

source("R/00_utils.R")
source("R/13_report.R")

out <- render_pdac_report(quiet = FALSE)
message("Report written to: ", out)
