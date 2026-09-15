# Data sources and provenance

## Discovery cohorts

| GEO accession | Primary role | Curated design |
|---|---|---|
| GSE15471 | Discovery | Paired PDAC tumor/normal; technical replicates averaged before modelling |
| GSE28735 | Discovery | Matched PDAC tumor/adjacent non-tumor pairs |
| GSE62165 | Discovery | PDAC and control |
| GSE16515 | Discovery | PDAC and normal; partly paired |
| GSE71989 | Discovery | PDAC, normal/control, chronic-pancreatitis sample retained only for metadata audit |
| GSE91035 | Transition analysis | Normal, benign, PDAC; two confirmed chronic-pancreatitis samples excluded from primary analyses |

The repository obtains GEO content with `GEOquery` and records expected-versus-observed metadata before modelling.

## TCGA / GTEx

Cross-platform validation uses the UCSC Xena Toil harmonized TCGA/TARGET/GTEx expression matrix and phenotype table. URLs and local cache paths are specified in `config/atlas_config.yml`.

The repository uses memory-safe gene/sample subsetting scripts so the full expression matrix does not need to be loaded into memory.

## Pan-cancer specificity

The specificity panel compares TCGA tumors with tissue-matched GTEx normals for:

- COAD
- STAD
- LIHC
- LUAD

Cancer identity is taken from disease metadata, not inferred from organ site alone.

## Single-cell reference

The public code expects a PDAC single-cell H5AD/R object or preconverted expression/metadata tables at the paths configured under `validation.scrna` in `config/atlas_config.yml`.

The full source object is not redistributed in this repository. Frozen single-cell-derived tables used in the release are included among the supplementary outputs.

## CIBERSORTx

The final composition validation used a seven-class PDAC single-cell reference:

- malignant epithelial
- ductal
- acinar
- fibroblast/stromal
- immune
- endothelial
- endocrine

CIBERSORTx itself is an external service/software product and is not redistributed here. The repository provides the reference/mixture export and result-import scripts; final frozen result tables are included under `results/frozen_release/`.

## DepMap and Open Targets

Optional target-annotation layers use local DepMap expression/CRISPR files and the Open Targets API. These resources are not redistributed. Relevant local file paths and switches are defined in `config/atlas_config.yml`.

## Frozen outputs

All supplementary tables and final validation summaries used for release v1.0.0 are stored under `results/frozen_release/`. These should be treated as the immutable paper-facing result set for this release.
