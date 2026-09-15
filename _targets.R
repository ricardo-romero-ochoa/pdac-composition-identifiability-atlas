library(targets)
library(tarchetypes)

source("R/00_utils.R")
source("R/01_download_geo.R")
source("R/02_metadata_curation.R")
source("R/03_preprocess_microarray.R")
source("R/04_de_limma.R")
source("R/05_meta_analysis.R")
source("R/06_concordance.R")
source("R/07_deconvolution_adjustment.R")
source("R/08_pathway_tf_activity.R")
source("R/09_transition_gse91035.R")
source("R/10_consensus_modules.R")
source("R/11_hubs_druggability.R")
source("R/12_figures.R")
source("R/13_report.R")
source("R/14_reinforcements.R")

tar_option_set(
  packages = c(
    "tidyverse", "data.table", "yaml", "janitor", "GEOquery", "Biobase",
    "limma", "metafor", "AnnotationDbi", "org.Hs.eg.db", "msigdbr",
    "GSVA", "pheatmap", "ggplot2", "ggrepel", "cowplot", "patchwork",
    "RColorBrewer", "matrixStats", "uwot", "WGCNA", "httr2", "jsonlite",
    "scales"
  ),
  format = "rds"
)

list(
  tar_target(cfg, read_config("config/atlas_config.yml")),
  tar_target(esets, download_all_geo(cfg), cue = tar_cue(mode = "thorough")),
  tar_target(metadata_list, curate_all_metadata(esets, cfg)),
  tar_target(metadata_audit, metadata_audit_report(metadata_list, cfg)),
  tar_target(dataset_list, preprocess_all(esets, metadata_list, cfg)),
  tar_target(dataset_summary, write_expression_summaries(dataset_list, cfg)),

  tar_target(de_tumor_vs_control, run_all_de(dataset_list, cfg, contrast = "tumor_vs_control")),
  tar_target(meta_gene, meta_analyze_de(de_tumor_vs_control, cfg, feature_col = "gene", out_prefix = "gene")),
  tar_target(core_signature, define_core_signature(meta_gene, cfg)),
  tar_target(tiered_signature, classify_signature_tiers(meta_gene, cfg)),
  tar_target(concordance, run_concordance(de_tumor_vs_control, cfg)),

  tar_target(deconv_scores, run_all_deconvolution(dataset_list, cfg)),
  tar_target(de_tme_adjusted, run_stroma_adjusted_de(dataset_list, deconv_scores, cfg)),
  tar_target(meta_gene_tme_adjusted, meta_analyze_de(de_tme_adjusted, cfg, feature_col = "gene", out_prefix = "gene_tme_adjusted")),
  tar_target(tme_adjustment_grid, run_tme_adjustment_grid(dataset_list, deconv_scores, cfg)),
  tar_target(gene_deconv_cor, run_gene_deconv_correlations(dataset_list, deconv_scores, cfg)),
  tar_target(microenvironment_signature, define_microenvironment_program(meta_gene, meta_gene_tme_adjusted, gene_deconv_cor, cfg)),
  tar_target(tumor_tme_classes, tme_retention_classification(meta_gene, meta_gene_tme_adjusted, gene_deconv_cor, cfg)),

  tar_target(activity_scores, run_all_activity(dataset_list, cfg)),
  tar_target(activity_results, run_all_activity_de_and_meta(activity_scores, dataset_list, cfg)),

  tar_target(transition_res, run_gse91035_transition(dataset_list, cfg)),
  tar_target(transition_signature, define_transition_program(transition_res, cfg)),
  tar_target(transition_validation, validate_transition_program_external(transition_signature, dataset_list, cfg)),

  tar_target(wgcna_res, run_consensus_wgcna(dataset_list, cfg)),
  tar_target(module_meta, meta_analyze_modules(wgcna_res, cfg)),
  tar_target(module_annotations, annotate_modules(wgcna_res, core_signature, microenvironment_signature, transition_signature, cfg)),

  tar_target(program_tables, write_program_tables(core_signature, microenvironment_signature, transition_signature, cfg)),
  tar_target(hub_table, make_integrative_hub_table(meta_gene, meta_gene_tme_adjusted, core_signature, microenvironment_signature, transition_signature, wgcna_res, module_meta, cfg)),
  tar_target(strict_hub_table, make_strict_integrative_hub_table(hub_table, tiered_signature, tumor_tme_classes, transition_validation, cfg)),
  tar_target(reinforcement_tables, write_reinforcement_tables(metadata_audit, tiered_signature, tumor_tme_classes, transition_validation, strict_hub_table, cfg)),

  tar_target(figures, make_all_figures(dataset_list, meta_gene, deconv_scores, transition_res, cfg)),
  tar_target(reinforcement_figures, make_reinforcement_figures(tiered_signature, tumor_tme_classes, transition_validation, strict_hub_table, cfg)),
  tar_target(report, {
    figures
    reinforcement_figures
    hub_table
    strict_hub_table
    program_tables
    reinforcement_tables
    render_pdac_report(quiet = TRUE)
  })
)
