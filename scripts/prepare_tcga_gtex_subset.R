#!/usr/bin/env Rscript
source("R/00_utils.R")
source("R/15_external_validation_utils.R")
source("R/16_tcga_gtex_validation.R")

cfg <- read_config("config/atlas_config.yml")
vcfg <- get_tcga_gtex_cfg(cfg)

cat("Step 1/2: exporting atlas signature gene universe...\n")
out_gene_file <- vcfg$requested_gene_file %||% "data/external/cache/tcga_gtex_requested_genes.txt"
dir.create(dirname(out_gene_file), recursive = TRUE, showWarnings = FALSE)
catalog <- build_signature_catalog(cfg)
genes <- signature_gene_universe(catalog)
genes <- sort(unique(genes[!is.na(genes) & nzchar(genes)]))
writeLines(genes, out_gene_file, useBytes = TRUE)
cat("Wrote ", length(genes), " genes to ", out_gene_file, "\n", sep = "")

cat("Step 2/2: streaming the UCSC Xena matrix with Python...\n")
python <- Sys.which("python3")
if (!nzchar(python)) python <- Sys.which("python")
if (!nzchar(python)) {
  stop("Neither python3 nor python was found on PATH. Install Python 3 or provide preprocessed TCGA/GTEx expression + metadata files.", call. = FALSE)
}

args <- c(
  "scripts/create_tcga_gtex_xena_subset.py",
  "--expression", vcfg$xena_expression_file %||% "data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz",
  "--phenotype", vcfg$xena_phenotype_file %||% "data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz",
  "--gene-list", out_gene_file,
  "--probemap", vcfg$xena_probemap_file %||% "data/external/tcga_gtex/gencode.v23.annotation.gene.probemap",
  "--probemap-url", vcfg$xena_probemap_url %||% "https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap",
  "--out-expression", vcfg$expression_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz",
  "--out-metadata", vcfg$metadata_file %||% "data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv",
  "--max-samples", as.character(vcfg$max_xena_samples %||% 1000L)
)
status <- system2(python, args = args)
if (!identical(status, 0L)) {
  stop("Python TCGA/GTEx subsetting failed with exit status ", status, call. = FALSE)
}
cat("TCGA/GTEx subset prepared. Now run:\n")
cat("  Rscript scripts/run_tcga_gtex_validation.R\n")
