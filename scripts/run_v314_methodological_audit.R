#!/usr/bin/env Rscript
# Runs the v3.14 methodological correction layer before report drafting.
# This intentionally rebuilds the base pipeline because ordinary limma SEs,
# within-tumor gene/deconvolution correlations, S1/S4 adjustment semantics, and
# transition identifier auditing change upstream tables.
if (!requireNamespace("targets", quietly = TRUE)) stop("Install the targets package first.", call. = FALSE)
targets::tar_make()
if (file.exists("_targets_external.R")) targets::tar_make(script = "_targets_external.R")
dir.create("results/session", recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "results/session/sessionInfo_v3_14.txt")
if (file.exists("scripts/inspect_v314_outputs.R")) source("scripts/inspect_v314_outputs.R")
