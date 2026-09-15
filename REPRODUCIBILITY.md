# Reproducing release v1.0.0

The repository contains both executable analysis code and a frozen copy of the paper-facing tables/figures.

## Core analysis

```bash
Rscript scripts/run_pipeline.R
Rscript scripts/run_v314_methodological_audit.R
Rscript scripts/run_reinforcements.R
```

## External validation

TCGA/GTEx:

```bash
Rscript scripts/export_tcga_gtex_gene_universe.R
Rscript scripts/prepare_tcga_gtex_subset.R
Rscript scripts/run_tcga_gtex_validation.R
```

Single-cell mapping:

```bash
Rscript scripts/prepare_scrna_tables.R
Rscript scripts/run_scrna_mapping.R
```

Pan-cancer and PurIST:

```bash
Rscript scripts/prepare_pancancer_xena_subset.R
Rscript scripts/run_pan_cancer_specificity_panel.R
Rscript scripts/run_purist_subtype_benchmark.R
```

Final CIBERSORTx validation:

```bash
Rscript scripts/export_cibersortx_final_validation_inputs.R
# run CIBERSORTx externally
Rscript scripts/import_cibersortx_final_validation.R
```

## Frozen outputs

The exact release tables and figures are already provided under:

```text
results/frozen_release/
```

This distinction is intentional: public upstream resources and web/API services can change, whereas the frozen release should remain auditable.

## Determinism

- analysis seed: `20260531`
- all primary thresholds: `config/atlas_config.yml`
- explicit sample overrides: `config/manual_overrides/sample_metadata_overrides.csv`
- release version: `1.0.0`
