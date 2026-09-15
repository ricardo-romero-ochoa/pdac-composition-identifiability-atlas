# Final reference-based composition validation with CIBERSORTx

## Why this is the final composition check

The definitive reference-based composition check uses CIBERSORTx with a PDAC single-cell reference. The public release reports this analysis as the final reference-based validation.

## Reference design

Seven broad PDAC compartments are prespecified:

1. `malignant_epithelial`
2. `ductal`
3. `acinar`
4. `fibroblast_stromal`
5. `immune`
6. `endothelial`
7. `endocrine`

Schwann/neural and unknown cells are excluded.

The exporter uses the full `scrna_reference.h5ad`, not the atlas-restricted
compact TSV. It prefers a counts layer or `adata.raw`, converts to CPM per cell,
uses at most 700 cells per class (<=4900 total), and retains up to 18,000
expressed genes.

## Step 1 — prepare reference and six mixture files

From the repository root:

```bash
Rscript scripts/export_cibersortx_final_validation_inputs.R
```

If Python lacks the h5ad dependencies:

```bash
python3 -m pip install --user anndata h5py scipy pandas numpy
```

The export directory is:

```text
results/external_deconvolution_inputs/cibersortx/
```

It contains:

- `PDAC_7class_refsample.txt`
- `PDAC_7class_reference_manifest.json`
- `PDAC_7class_reference_cells.tsv`
- `PDAC_7class_reference_class_counts.tsv`
- six files under `mixtures/`

The mixture exporter uses the cached, curated six-cohort matrices and averages
technical replicates consistently with the atlas pipeline. Log2-scale microarray
matrices are returned to linear scale with `2^x`.

## Step 2 — run CIBERSORTx

CIBERSORTx requires an academic/non-commercial account and token. Never commit
the token.

In WSL/Linux:

```bash
export CIBERSORTX_EMAIL='your_academic_email'
export CIBERSORTX_TOKEN='your_token'
bash scripts/run_cibersortx_final_validation.sh
```

The runner uses the official `cibersortx/fractions` Docker image and runs each
cohort independently with:

- single-cell reference mode;
- S-mode batch correction;
- quantile normalization enabled for these microarray mixtures;
- relative fractions (`absolute=FALSE`);
- 50 permutations.

Each cohort result is copied to:

```text
results/external_deconvolution_inputs/cibersortx/results/
```

### Web-interface alternative

If Docker is unavailable, use the same reference and each mixture file in the
CIBERSORTx fractions workflow. Use single-cell reference mode, S-mode batch
correction, relative fractions, and quantile normalization for the microarray
mixtures. Download each fractions result and name it:

```text
GSE15471_CIBERSORTx_Adjusted.txt
GSE28735_CIBERSORTx_Adjusted.txt
GSE62165_CIBERSORTx_Adjusted.txt
GSE16515_CIBERSORTx_Adjusted.txt
GSE71989_CIBERSORTx_Adjusted.txt
GSE91035_CIBERSORTx_Adjusted.txt
```

Place them in:

```text
results/external_deconvolution_inputs/cibersortx/results/
```

Then run:

```bash
Rscript scripts/import_cibersortx_final_validation.R
```

## Primary validation comparisons

The analysis is intentionally within-cohort first, avoiding pooled
cross-platform offsets:

- `fibroblast_stromal` vs `stromal_caf`
- `immune` vs `immune_pan`
- `endothelial` vs `endothelial`
- `acinar` vs `acinar_pancreas`
- `ductal` vs `ductal_epithelial`

Outputs:

- `reference_deconvolution_expected_pair_by_dataset.tsv`
- `reference_deconvolution_expected_pair_meta.tsv`
- `reference_deconvolution_qc.tsv`
- `reference_deconvolution_summary.tsv`
- `reference_deconvolution_marker_score_joined.tsv`

The current readiness policy accepts only a summary containing a
`CIBERSORTx` method. The exploratory local NNLS result can no longer make the
reference-deconvolution gate appear complete.

## Scientific decision rule

Do not invent a new post-hoc numeric threshold merely to obtain a PASS.
Evaluate:

1. direction and consistency of the five expected biological pairings;
2. within-study Spearman correlations across the six cohorts;
3. CIBERSORTx fit diagnostics;
4. whether immune/endothelial concordance improves relative to the local
   rank-simplex sensitivity analysis.

If CIBERSORTx still contradicts immune/endothelial marker scores, retain the
result and weaken the corresponding composition claims rather than tuning
another deconvolution model.

The archive DOI remains the last administrative action and should be minted
only after this result and the manuscript revision are frozen.
