.PHONY: install base audit reinforce tcga-gtex scrna pancancer purist cibersortx-inputs reports

install:
	Rscript scripts/install_packages.R

base:
	Rscript scripts/run_pipeline.R

audit:
	Rscript scripts/run_v314_methodological_audit.R

reinforce:
	Rscript scripts/run_reinforcements.R

tcga-gtex:
	Rscript scripts/export_tcga_gtex_gene_universe.R
	Rscript scripts/prepare_tcga_gtex_subset.R
	Rscript scripts/run_tcga_gtex_validation.R

scrna:
	Rscript scripts/prepare_scrna_tables.R
	Rscript scripts/run_scrna_mapping.R

pancancer:
	Rscript scripts/prepare_pancancer_xena_subset.R
	Rscript scripts/run_pan_cancer_specificity_panel.R

purist:
	Rscript scripts/run_purist_subtype_benchmark.R

cibersortx-inputs:
	Rscript scripts/export_cibersortx_final_validation_inputs.R

reports:
	Rscript scripts/render_report.R
	Rscript scripts/render_external_validation_report.R
