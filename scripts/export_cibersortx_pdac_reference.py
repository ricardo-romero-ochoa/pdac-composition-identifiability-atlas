#!/usr/bin/env python3
"""
Export a 7-class PDAC single-cell reference in CIBERSORTx refsample format.

Primary classes:
  malignant_epithelial
  ductal
  acinar
  fibroblast_stromal
  immune
  endothelial
  endocrine

The exporter deliberately uses the full h5ad reference rather than the compact
atlas-gene TSV. It prefers a true counts layer, then adata.raw, then count-like
adata.X. If only log1p-normalized X is available and adata.uns['log1p'] is
present, it reconstructs linear values with expm1 and records that choice.
Scaled/negative matrices are rejected.

Output refsample format follows the CIBERSORTx single-cell convention:
the first row contains repeated cell-class labels; subsequent rows contain
GeneSymbol and expression values. No header row is added above the class row.
"""

from __future__ import annotations
import argparse
import csv
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

try:
    import anndata as ad
    import scipy.sparse as sp
except Exception as exc:
    raise SystemExit(
        "Python dependencies are missing. Install them in the Python environment "
        "used for this script:\n"
        "  python3 -m pip install --user anndata h5py scipy pandas numpy\n"
        f"Original import error: {exc}"
    )

UNKNOWN = {"", "unknown", "unk", "na", "n/a", "nan", "none", "null", "unassigned"}

ANNOTATION_PRIORITY = [
    "new_celltypes",
    "celltype_infercnv",
    "Level 1 Annotation",
    "broad_celltypes",
    "cell_type",
    "celltype",
    "annotation",
    "cell_annotation",
    "final_annotation",
]

GENE_SYMBOL_COLUMNS = [
    "gene_symbol", "gene_symbols", "symbol", "hgnc_symbol",
    "GeneSymbol", "gene_name", "gene_names",
]

COUNTS_LAYERS = ["counts", "raw_counts", "count", "Counts", "raw"]

TARGET_CLASSES = [
    "malignant_epithelial",
    "ductal",
    "acinar",
    "fibroblast_stromal",
    "immune",
    "endothelial",
    "endocrine",
]


def norm_text(x: object) -> str:
    return re.sub(r"\s+", " ", str(x).strip().lower())


def map_pdac7(label: object) -> str | None:
    z = norm_text(label)
    if z in UNKNOWN:
        return None

    # Specific epithelial subclasses first.
    if re.search(r"malignan|cancer|tumou?r", z):
        return "malignant_epithelial"
    if re.search(r"ductal|atypical[_ -]?duct|pancreatic duct", z):
        return "ductal"
    if re.search(r"acinar", z):
        return "acinar"

    if re.search(r"fibro|caf|stellate|pericyte|smooth muscle|vsmc|strom|mesench", z):
        return "fibroblast_stromal"
    if re.search(r"endothel|capillary|arterial|venous|vascular", z):
        return "endothelial"
    if re.search(
        r"immune|myeloid|macroph|monocyte|lymph|t[ _-]?cell|b[ _-]?cell|"
        r"\bnk\b|mast|dendritic|neutroph|plasma", z
    ):
        return "immune"
    if re.search(r"endocrine|\bbeta\b|\balpha\b|\bdelta\b|\bgamma\b|islet", z):
        return "endocrine"
    return None


def choose_annotation(obs) -> tuple[str, list[str | None], list[dict]]:
    cols = list(obs.columns)
    lower = {c.lower(): c for c in cols}
    candidates = []
    for c in ANNOTATION_PRIORITY:
        if c.lower() in lower:
            candidates.append(lower[c.lower()])
    for c in cols:
        lc = c.lower()
        if c not in candidates and re.search(r"cell.*type|type.*cell|annot|compartment|lineage", lc):
            candidates.append(c)

    if not candidates:
        raise SystemExit(
            "No plausible cell-type annotation column found in h5ad obs. "
            f"Available columns: {', '.join(cols)}"
        )

    diagnostics = []
    best = None
    best_key = None
    best_mapped = None
    for order, col in enumerate(candidates):
        mapped = [map_pdac7(v) for v in obs[col].tolist()]
        counts = Counter(x for x in mapped if x is not None)
        n_ge20 = sum(v >= 20 for v in counts.values())
        n_target = sum(1 for cls in TARGET_CLASSES if counts.get(cls, 0) >= 20)
        n_mapped = sum(counts.values())
        key = (n_target, n_ge20, n_mapped, -order)
        diagnostics.append({
            "annotation_column": col,
            "candidate_order": order + 1,
            "n_mapped": n_mapped,
            "n_target_classes_ge20": n_target,
            "mapped_counts": ";".join(f"{k}={counts.get(k,0)}" for k in TARGET_CLASSES),
        })
        if best_key is None or key > best_key:
            best = col
            best_key = key
            best_mapped = mapped

    assert best is not None and best_mapped is not None
    counts = Counter(x for x in best_mapped if x is not None)
    missing = [cls for cls in TARGET_CLASSES if counts.get(cls, 0) < 20]
    if missing:
        raise SystemExit(
            "No annotation column yields all seven PDAC reference classes with >=20 cells. "
            f"Best column: {best}; counts: "
            + ", ".join(f"{k}={counts.get(k,0)}" for k in TARGET_CLASSES)
            + ". Missing/underpowered: " + ", ".join(missing)
        )
    return best, best_mapped, diagnostics


def choose_gene_symbols(var, var_names) -> np.ndarray:
    lower = {str(c).lower(): c for c in var.columns}
    for c in GENE_SYMBOL_COLUMNS:
        if c.lower() in lower:
            vals = var[lower[c.lower()]].astype(str).to_numpy()
            if np.mean([bool(v and v.lower() not in UNKNOWN) for v in vals]) > 0.8:
                return vals
    return np.asarray([str(x) for x in var_names], dtype=object)


def looks_count_like(X) -> bool:
    # Sample a bounded number of finite values.
    if sp.issparse(X):
        vals = np.asarray(X.data)
    else:
        vals = np.asarray(X).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return False
    if vals.size > 200000:
        rng = np.random.default_rng(314150)
        vals = rng.choice(vals, 200000, replace=False)
    if np.nanmin(vals) < -1e-8:
        return False
    frac_integerish = np.mean(np.abs(vals - np.round(vals)) < 1e-6)
    q99 = float(np.quantile(vals, 0.99))
    return frac_integerish > 0.80 or q99 > 25


def choose_expression_source(adata):
    for layer in COUNTS_LAYERS:
        if layer in adata.layers:
            return adata.layers[layer], adata.var, adata.var_names, f"layer:{layer}", False

    if adata.raw is not None:
        return adata.raw.X, adata.raw.var, adata.raw.var_names, "adata.raw.X", False

    if looks_count_like(adata.X):
        return adata.X, adata.var, adata.var_names, "adata.X_count_like", False

    # Explicit log1p metadata is enough to permit a reversible reconstruction.
    if "log1p" in adata.uns:
        X = adata.X
        if sp.issparse(X):
            if X.data.size and np.nanmin(X.data) < -1e-8:
                raise SystemExit(
                    "adata.X contains negative values and appears scaled; cannot reconstruct "
                    "a linear CIBERSORTx reference. Provide a counts layer or adata.raw."
                )
        else:
            if np.nanmin(np.asarray(X)) < -1e-8:
                raise SystemExit(
                    "adata.X contains negative values and appears scaled; cannot reconstruct "
                    "a linear CIBERSORTx reference. Provide a counts layer or adata.raw."
                )
        return X, adata.var, adata.var_names, "adata.X_log1p_reconstructed", True

    raise SystemExit(
        "Could not identify a non-log/count-like expression source in the h5ad. "
        "Expected a counts/raw_counts layer, adata.raw, count-like adata.X, or "
        "adata.uns['log1p']. No scaled/negative matrix will be exported."
    )


def balanced_cell_selection(mapped, max_per_class: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    keep = []
    mapped_arr = np.asarray(mapped, dtype=object)
    for cls in TARGET_CLASSES:
        idx = np.where(mapped_arr == cls)[0]
        if idx.size == 0:
            continue
        if idx.size > max_per_class:
            idx = rng.choice(idx, max_per_class, replace=False)
        keep.extend(idx.tolist())
    return np.asarray(sorted(keep), dtype=int)


def subset_gene_cell_matrix(X, gene_idx, cell_idx, n_obs, n_vars):
    """
    Return a genes x cells matrix from an AnnData expression source.

    Standard AnnData convention is cells x genes (n_obs x n_vars). Downstream
    CIBERSORTx export code in this script works in genes x cells orientation,
    so slice cells first, genes second, then transpose.

    A defensive already-transposed branch is retained for unusual inputs.
    """
    shape = tuple(X.shape)

    if shape == (n_obs, n_vars):
        if sp.issparse(X):
            return X[cell_idx, :][:, gene_idx].T.tocsr()
        return np.asarray(X)[np.ix_(cell_idx, gene_idx)].T

    if shape == (n_vars, n_obs):
        if sp.issparse(X):
            return X[gene_idx, :][:, cell_idx].tocsr()
        return np.asarray(X)[np.ix_(gene_idx, cell_idx)]

    raise ValueError(
        "Unexpected expression matrix shape "
        f"{shape}; expected AnnData cells x genes {(n_obs, n_vars)} "
        f"or transposed genes x cells {(n_vars, n_obs)}."
    )


def column_sums(X):
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel()
    return np.asarray(X).sum(axis=0)


def row_nonzero_fraction(X):
    if sp.issparse(X):
        return np.asarray((X > 0).sum(axis=1)).ravel() / X.shape[1]
    return np.mean(np.asarray(X) > 0, axis=1)


def row_means(X):
    if sp.issparse(X):
        return np.asarray(X.mean(axis=1)).ravel()
    return np.asarray(X).mean(axis=1)


def to_dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def sanitize_symbol(s: object) -> str:
    x = str(s).strip()
    if "|" in x:
        parts = [p for p in x.split("|") if p]
        # Prefer a symbol-like non-Ensembl token.
        non_ens = [p for p in parts if not re.match(r"^ENSG\d+", p, flags=re.I)]
        if non_ens:
            x = non_ens[-1]
    x = re.sub(r"\.\d+$", "", x) if re.match(r"^ENSG\d+\.\d+$", x, flags=re.I) else x
    return x.upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/external/scrna/scrna_reference.h5ad")
    ap.add_argument("--outdir", default="results/external_deconvolution_inputs/cibersortx")
    ap.add_argument("--max-cells-per-class", type=int, default=700)
    ap.add_argument("--min-expression-fraction", type=float, default=0.01)
    ap.add_argument("--max-genes", type=int, default=18000)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()

    h5ad = Path(args.h5ad)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not h5ad.exists():
        raise SystemExit(f"Full scRNA h5ad not found: {h5ad}")

    print(f"Reading full scRNA reference: {h5ad}")
    adata = ad.read_h5ad(h5ad)

    ann_col, mapped, ann_diag = choose_annotation(adata.obs)
    selected_cells = balanced_cell_selection(mapped, args.max_cells_per_class, args.seed)
    selected_labels = np.asarray(mapped, dtype=object)[selected_cells]
    selected_names = np.asarray(adata.obs_names.astype(str))[selected_cells]

    X0, var, var_names, source_name, reconstruct_log1p = choose_expression_source(adata)
    symbols = choose_gene_symbols(var, var_names)
    symbols = np.asarray([sanitize_symbol(x) for x in symbols], dtype=object)

    # Restrict to selected cells first.
    gene_idx = np.arange(len(symbols), dtype=int)
    X = subset_gene_cell_matrix(
        X0,
        gene_idx=gene_idx,
        cell_idx=selected_cells,
        n_obs=adata.n_obs,
        n_vars=len(var_names),
    )

    if X.shape != (len(gene_idx), len(selected_cells)):
        raise RuntimeError(
            "Internal orientation invariant failed after AnnData slicing: "
            f"got {X.shape}, expected {(len(gene_idx), len(selected_cells))} "
            "(genes x cells)."
        )

    if reconstruct_log1p:
        if sp.issparse(X):
            X = X.copy()
            X.data = np.expm1(X.data)
        else:
            X = np.expm1(np.asarray(X))

    # Keep usable symbols only.
    valid_symbol = np.asarray([
        bool(s) and s.lower() not in UNKNOWN and not s.startswith("ENSG")
        for s in symbols
    ])
    gene_idx = np.where(valid_symbol)[0]
    X = X[gene_idx, :] if sp.issparse(X) else X[gene_idx, :]
    symbols = symbols[gene_idx]

    # Deduplicate symbols by keeping the row with highest mean expression.
    means = row_means(X)
    best = {}
    for i, (sym, mu) in enumerate(zip(symbols, means)):
        if sym not in best or mu > best[sym][1]:
            best[sym] = (i, float(mu))
    keep = np.asarray(sorted(v[0] for v in best.values()), dtype=int)
    X = X[keep, :] if sp.issparse(X) else X[keep, :]
    symbols = symbols[keep]

    # Normalize cells to CPM after reconstruction/count selection.
    libs = column_sums(X)
    good_cells = np.isfinite(libs) & (libs > 0)
    if not np.all(good_cells):
        X = X[:, good_cells]
        selected_labels = selected_labels[good_cells]
        selected_names = selected_names[good_cells]
        libs = libs[good_cells]

    if sp.issparse(X):
        X = X @ sp.diags(1e6 / libs)
    else:
        X = np.asarray(X) * (1e6 / libs)[None, :]

    # Filter low-information genes and cap to keep the upload tractable.
    frac = row_nonzero_fraction(X)
    keep = np.where(frac >= args.min_expression_fraction)[0]
    X = X[keep, :] if sp.issparse(X) else X[keep, :]
    symbols = symbols[keep]

    means = row_means(X)
    if len(symbols) > args.max_genes:
        top = np.argsort(means)[-args.max_genes:]
        top = np.sort(top)
        X = X[top, :] if sp.issparse(X) else X[top, :]
        symbols = symbols[top]

    # Final ordering by gene symbol makes the file deterministic.
    ord_idx = np.argsort(symbols)
    X = X[ord_idx, :] if sp.issparse(X) else X[ord_idx, :]
    symbols = symbols[ord_idx]

    counts = Counter(selected_labels.tolist())
    if any(counts.get(cls, 0) < 20 for cls in TARGET_CLASSES):
        raise SystemExit(
            "The final exported reference unexpectedly lost a required class: "
            + ", ".join(f"{c}={counts.get(c,0)}" for c in TARGET_CLASSES)
        )

    ref_file = outdir / "PDAC_7class_refsample.txt"
    dense = to_dense(X).astype(np.float64, copy=False)

    print(
        f"Writing CIBERSORTx refsample: {len(symbols)} genes x "
        f"{dense.shape[1]} cells -> {ref_file}"
    )
    with ref_file.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        # CIBERSORTx refsample first row = class labels, not cell IDs.
        w.writerow(["GeneSymbol"] + selected_labels.tolist())
        for i, gene in enumerate(symbols):
            vals = [format(float(v), ".7g") for v in dense[i, :]]
            w.writerow([gene] + vals)

    # Metadata and diagnostics.
    meta_file = outdir / "PDAC_7class_reference_cells.tsv"
    obs_sel = adata.obs.iloc[selected_cells].copy()
    original_labels = obs_sel[ann_col].astype(str).tolist()
    donor_col = next(
        (c for c in ["sample_id", "sampleid", "pid", "donor_id", "patient_id", "donor", "patient"]
         if c in adata.obs.columns),
        None,
    )
    donors = obs_sel[donor_col].astype(str).tolist() if donor_col else ["NA"] * len(obs_sel)

    with meta_file.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cell", "cibersortx_class", "original_annotation", "annotation_column", "donor"])
        for cell, cls, lab, donor in zip(selected_names, selected_labels, original_labels, donors):
            w.writerow([cell, cls, lab, ann_col, donor])

    diag_file = outdir / "PDAC_7class_annotation_candidates.tsv"
    with diag_file.open("w", encoding="utf-8", newline="") as fh:
        fields = ["annotation_column", "candidate_order", "n_mapped", "n_target_classes_ge20", "mapped_counts"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(ann_diag)

    class_file = outdir / "PDAC_7class_reference_class_counts.tsv"
    with class_file.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cibersortx_class", "n_cells"])
        for cls in TARGET_CLASSES:
            w.writerow([cls, counts.get(cls, 0)])

    manifest = {
        "h5ad": str(h5ad),
        "annotation_column": ann_col,
        "expression_source": source_name,
        "log1p_reconstructed": bool(reconstruct_log1p),
        "normalization": "CPM_per_cell",
        "n_cells": int(dense.shape[1]),
        "n_genes": int(dense.shape[0]),
        "max_cells_per_class": int(args.max_cells_per_class),
        "min_expression_fraction": float(args.min_expression_fraction),
        "max_genes": int(args.max_genes),
        "seed": int(args.seed),
        "class_counts": {cls: int(counts.get(cls, 0)) for cls in TARGET_CLASSES},
    }
    (outdir / "PDAC_7class_reference_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print("Reference export complete.")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
