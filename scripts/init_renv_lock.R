#!/usr/bin/env Rscript
# Create a real renv.lock on the machine that runs the analysis. This script is
# intentionally not replaced by a fabricated lockfile because package versions
# must reflect the executing R/Bioconductor environment.
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init(bare = TRUE, restart = FALSE)
if (requireNamespace("BiocManager", quietly = TRUE)) {
  message("Bioconductor version: ", as.character(BiocManager::version()))
}
renv::snapshot(prompt = FALSE)
