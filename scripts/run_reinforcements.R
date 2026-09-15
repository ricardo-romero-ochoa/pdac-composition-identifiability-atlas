#!/usr/bin/env Rscript
# Run only the manuscript-reinforcement layer and its dependencies.
if (!requireNamespace("targets", quietly = TRUE)) stop("Install targets first: install.packages('targets')", call. = FALSE)
targets::tar_make(names = c(
  "metadata_audit",
  "tiered_signature",
  "tumor_tme_classes",
  "transition_validation",
  "strict_hub_table",
  "reinforcement_tables",
  "reinforcement_figures",
  "report"
))
