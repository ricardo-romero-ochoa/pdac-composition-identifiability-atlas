#!/usr/bin/env Rscript

cran <- c(
  "tidyverse", "data.table", "yaml", "scales", "janitor", "metafor", "pheatmap",
  "ggrepel", "cowplot", "patchwork", "RColorBrewer", "uwot", "matrixStats",
  "circlize", "httr2", "jsonlite", "rmarkdown", "knitr", "targets",
  "tarchetypes", "remotes", "survival", "survminer"
)

bioc <- c(
  "GEOquery", "Biobase", "limma", "AnnotationDbi", "org.Hs.eg.db",
  "msigdbr", "GSVA", "ComplexHeatmap", "WGCNA", "hgu133plus2.db",
  "hugene10sttranscriptcluster.db", "pd.hugene.1.0.st.v1", "progeny",
  "dorothea", "viper", "SummarizedExperiment", "SingleCellExperiment", "zellkonverter"
)

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) installer(missing)
}

install_if_missing(cran, function(x) install.packages(x, repos = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
install_if_missing(bioc, function(x) BiocManager::install(x, ask = FALSE, update = FALSE))

# Optional packages. Failures are non-fatal because some are not always available
# for all R/Bioconductor combinations.
optional <- c("xCell", "immunedeconv", "STRINGdb", "Seurat", "SeuratObject", "UCSCXenaTools", "TCGAbiolinks")
for (pkg in optional) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Optional package not installed: ", pkg, ". Install manually if needed.")
  }
}

message("Package installation check completed.")
