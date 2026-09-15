# Composition identifiability in public PDAC transcriptomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Release](https://img.shields.io/badge/release-v1.0.0-blue.svg)](CITATION.cff)

Reproducible code and frozen analysis outputs for:

**Composition identifiability in public PDAC transcriptomics: a reproducible cross-cohort re-analysis and annotation resource**

This repository analyzes six public pancreatic ductal adenocarcinoma (PDAC) transcriptomic cohorts and explicitly separates **cross-cohort reproducibility** from **cellular composition, subtype association, and organ specificity**. It is intended as a reproducible analysis/resource repository rather than a de novo biomarker-discovery package.

## Scientific scope

The discovery atlas uses:

- **GSE15471**
- **GSE28735**
- **GSE62165**
- **GSE16515**
- **GSE71989**
- **GSE91035**

The workflow includes metadata audit, paired or mixed-design differential expression, ordinary-standard-error random-effects meta-analysis, leave-one-dataset-out stability, a pre-specified S0-S4 composition-adjustment grid, within-group gene/composition correlations, transition-program analysis, pathway/regulator summaries, and external validation.

The frozen release also contains the final validation layers used in the associated study:

- **TCGA-PAAD / GTEx pancreas** module validation and TCGA survival analysis;
- **single-cell localization** of atlas programs;
- **pan-cancer specificity testing** against COAD, STAD, LIHC and LUAD with tissue-matched GTEx normals;
- **PurIST** sample-level basal-like/classical classification in 178 TCGA-PAAD tumors;
- **CIBERSORTx** 7-class reference-based deconvolution across 414 samples from the six discovery cohorts;
- prior-art overlap, subtype-signature overlap, DepMap and Open Targets annotation.

The main interpretation is deliberately conservative: reproducibility in bulk PDAC does **not** by itself establish tumor-cell autonomy or pancreas specificity.

## Frozen release results

Paper-facing tables and figures are versioned under:

```text
results/frozen_release/
├── supplementary_tables/
├── final_validation/
├── main_figures/
└── supplementary_figures/
```

These files correspond to release **v1.0.0** and are included so the archived repository remains directly auditable even when external services change.

Selected final checks include:

- pan-cancer testing showing strong transfer of the Tier-1 program to colorectal and gastric adenocarcinoma;
- PurIST association strongest for the transition composite (AUC 0.722), while the composition-adjusted retained program is subtype-neutral (absolute AUC 0.518);
- CIBERSORTx validation with median within-cohort Spearman correlations of **0.955** (acinar), **0.931** (stromal/CAF), **0.875** (immune) and **0.663** (endothelial).

## Repository structure

```text
R/                      Core analysis functions
scripts/                Reproducible entry points and data-preparation utilities
config/                 Analysis configuration and explicit metadata overrides
resources/              Curated literature/subtype/druggability resources
reports/                 R Markdown analysis reports
docs/                    Public reproducibility and validation documentation
results/frozen_release/  Frozen tables and figures for release v1.0.0
```

Raw GEO, Xena, DepMap and single-cell source files are intentionally not redistributed here. They are public or externally distributed resources and are downloaded/provided separately as documented below.

## Quick start

### 1. Software

R >= 4.3 is recommended. Install the R dependencies listed in `DESCRIPTION`:

```bash
Rscript scripts/install_packages.R
```

For H5AD processing and the CIBERSORTx reference exporter:

```bash
python3 -m pip install -r requirements-python.txt
```

A `Dockerfile` is provided for a container-oriented R environment.

### 2. Run the six-cohort GEO atlas

```bash
Rscript scripts/run_pipeline.R
```

The first run downloads the configured GEO expression/annotation resources and writes intermediate data under `data/` and analysis outputs under `results/`.

The two confirmed chronic-pancreatitis samples in GSE91035 are explicitly excluded through:

```text
config/manual_overrides/sample_metadata_overrides.csv
```

### 3. Run the methodological audit and reinforcement layer

```bash
Rscript scripts/run_v314_methodological_audit.R
Rscript scripts/run_reinforcements.R
```

### 4. TCGA/GTEx validation

The memory-safe UCSC Xena subset workflow is:

```bash
Rscript scripts/export_tcga_gtex_gene_universe.R
Rscript scripts/prepare_tcga_gtex_subset.R
Rscript scripts/run_tcga_gtex_validation.R
```

Large Xena source matrices are not stored in Git.

### 5. Pan-cancer specificity and PurIST validation

After the Xena source files are locally available:

```bash
Rscript scripts/prepare_pancancer_xena_subset.R
Rscript scripts/run_pan_cancer_specificity_panel.R
Rscript scripts/run_purist_subtype_benchmark.R
```

### 6. Single-cell mapping

Provide the PDAC single-cell reference as configured in `config/atlas_config.yml`, then:

```bash
Rscript scripts/prepare_scrna_tables.R
Rscript scripts/run_scrna_mapping.R
```

### 7. Final CIBERSORTx validation

Generate the seven-class PDAC single-cell reference and the six bulk mixture files:

```bash
Rscript scripts/export_cibersortx_final_validation_inputs.R
```

Run CIBERSORTx externally (web interface or the documented Docker workflow), then place the returned cohort results under:

```text
results/external_deconvolution_inputs/cibersortx/results/
```

and import them with:

```bash
Rscript scripts/import_cibersortx_final_validation.R
```

See [`docs/cibersortx_final_validation.md`](docs/cibersortx_final_validation.md).

## Reproducibility notes

- All primary thresholds and dataset-specific metadata rules are in `config/atlas_config.yml`.
- Confirmed manual exclusions are version controlled in `config/manual_overrides/`.
- The random-effects meta-analysis uses ordinary limma standard errors rather than moderated-t-derived pseudo-standard-errors.
- Composition adjustment is treated as a **specification grid**, not as proof of a uniquely correct composition-free model.
- `S1` (stromal + immune) is the primary parsimonious adjustment; epithelial/acinar adjustment is retained as an over-adjustment sensitivity analysis.
- The final CIBERSORTx analysis is the reference-based composition validation. Earlier exploratory local rank-simplex deconvolution code is intentionally not part of the public release.
- Internet-facing annotations can change over time; the frozen paper-facing outputs preserve the release used for the associated study.

## Data sources

See [`DATA_SOURCES.md`](DATA_SOURCES.md) for accession numbers, external-resource roles, and local file expectations.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). When this release is archived in Zenodo or another DOI-minting repository, please cite both the associated article and the archived software/data release.

## License

Code is released under the [MIT License](LICENSE). External datasets remain subject to the terms of their original repositories/providers.
