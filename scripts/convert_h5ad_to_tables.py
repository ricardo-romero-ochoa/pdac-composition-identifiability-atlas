#!/usr/bin/env python3
"""Memory-aware conversion of a PDAC scRNA/snRNA h5ad file to compact TSV tables.

The R targets step can then read the TSV files without requiring anndata inside
reticulate/zellkonverter. The converter keeps only atlas-relevant genes and can
optionally downsample cells per cell type.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd

try:
    import anndata as ad
except Exception as e:  # pragma: no cover
    sys.stderr.write(
        "Could not import anndata in this Python interpreter: %s\n"
        "Install into this exact Python with:\n  %s -m pip install anndata h5py pandas scipy numpy\n"
        % (e, sys.executable)
    )
    raise

try:
    import scipy.sparse as sp
except Exception:  # pragma: no cover
    sp = None


def norm_gene(x: object) -> str:
    s = "" if x is None else str(x)
    s = re.sub(r"^.*\\|", "", s)
    s = re.sub(r"\\..*$", "", s)
    return s.strip().upper()


def read_gene_universe(path: str) -> List[str]:
    if not path or not os.path.exists(path):
        return []
    tab = pd.read_csv(path, sep="\t", dtype=str)
    col = "gene" if "gene" in tab.columns else tab.columns[0]
    genes = [norm_gene(x) for x in tab[col].dropna().tolist()]
    return sorted(set(g for g in genes if g))


def infer_first_existing(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def make_symbol_index(adata) -> pd.DataFrame:
    var = adata.var.copy()
    var_names = pd.Series(list(adata.var_names), index=var.index, dtype="object")
    sym_col = infer_first_existing(
        var.columns,
        [
            "gene_symbols",
            "gene_symbol",
            "symbol",
            "feature_name",
            "gene_name",
            "name",
            "external_gene_name",
        ],
    )
    symbols = var[sym_col].astype(str) if sym_col is not None else var_names.astype(str)
    out = pd.DataFrame(
        {
            "var_pos": np.arange(adata.n_vars, dtype=int),
            "var_name": list(adata.var_names),
            "gene": [norm_gene(x) for x in symbols.tolist()],
        }
    )
    out = out[out["gene"].astype(bool)].copy()
    out = out.drop_duplicates("gene", keep="first")
    return out


UNKNOWN_TOKENS = {"", "unknown", "unk", "na", "n/a", "nan", "none", "null", "unspecified", "unassigned", "unannotated"}


def useful_column(meta: pd.DataFrame, col: str) -> bool:
    vals = meta[col].astype(str).str.strip().str.lower()
    vals = vals[~vals.isin(UNKNOWN_TOKENS)]
    return vals.nunique(dropna=True) > 0


def infer_metadata_column(meta: pd.DataFrame, explicit: Optional[str], candidates: List[str]) -> Optional[str]:
    if explicit:
        hit = infer_first_existing(meta.columns, [explicit])
        if hit is None:
            raise SystemExit(
                "Configured scRNA metadata column was not found: %s\nAvailable columns: %s"
                % (explicit, ", ".join(map(str, meta.columns)))
            )
        return hit
    for cand in candidates:
        hit = infer_first_existing(meta.columns, [cand])
        if hit is not None and useful_column(meta, hit):
            return hit
    bad = re.compile(r"barcode|cell$|sample|patient|donor|library|batch|ncount|nfeature|percent|mito|doublet|sex|age|stage|grade|score|leiden|cluster")
    good = re.compile(r"cell|type|annot|label|class|compartment|lineage|subtype|identity")
    for col in meta.columns:
        low = str(col).lower()
        if good.search(low) and not bad.search(low) and useful_column(meta, col):
            n = meta[col].astype(str).str.strip().nunique(dropna=True)
            if 2 <= n <= 100:
                return col
    return None


def infer_metadata(
    obs: pd.DataFrame,
    cell_ids: List[str],
    cell_type_column: Optional[str] = None,
    condition_column: Optional[str] = None,
    sample_column: Optional[str] = None,
    fail_if_unknown_celltype: bool = True,
    out_dir: Optional[str] = None,
) -> pd.DataFrame:
    meta = obs.copy()
    meta.insert(0, "cell", cell_ids)
    ct_candidates = [
        "cell_type", "cell_types", "celltype", "cell type", "CellType", "Cell_type", "cellType",
        "celltype_major", "cell_type_major", "major_cell_type", "major_celltype",
        "broad_cell_type", "broad_celltype", "cell_type_broad", "annotation", "annotations",
        "cell_annotation", "manual_annotation", "final_annotation", "author_cell_type",
        "author_celltype", "predicted_cell_type", "predicted.celltype", "predicted.id",
        "predicted_ID", "cell_ontology_class", "cell_ontology", "CellOntology", "subclass",
        "class", "compartment", "lineage", "cell.labels", "cell_label", "labels", "ident",
        "seurat_clusters",
    ]
    cond_candidates = [
        "condition", "condition_original", "disease", "Disease", "diagnosis", "group", "tissue",
        "sample_type", "tumor_normal", "status", "phenotype", "pathology",
    ]
    sample_candidates = [
        "sample", "sample_id", "Sample", "SampleID", "patient", "patient_id", "donor", "donor_id",
        "orig.ident", "library", "library_id", "specimen", "case", "case_id",
    ]

    ct = infer_metadata_column(meta, cell_type_column, ct_candidates)
    cond = infer_metadata_column(meta, condition_column, cond_candidates)
    sid = infer_metadata_column(meta, sample_column, sample_candidates)

    meta["cell_type"] = meta[ct].astype(str).str.strip() if ct is not None else "unknown"
    meta.loc[meta["cell_type"].str.lower().isin(UNKNOWN_TOKENS), "cell_type"] = "unknown"
    meta["condition"] = meta[cond].astype(str).str.strip() if cond is not None else "unknown"
    meta["sample_id"] = meta[sid].astype(str).str.strip() if sid is not None else meta["cell"]

    if out_dir:
        diag = []
        for col in meta.columns:
            vals = meta[col].astype(str).str.strip()
            examples = "; ".join(vals[vals != ""].drop_duplicates().head(8).tolist())
            diag.append({
                "column": col,
                "n_non_missing": int((vals != "").sum()),
                "n_unique": int(vals.nunique(dropna=True)),
                "example_values": examples,
            })
        pd.DataFrame(diag).to_csv(os.path.join(out_dir, "scrna_metadata_columns.tsv"), sep="\t", index=False)

    if fail_if_unknown_celltype and (meta["cell_type"] != "unknown").sum() == 0:
        raise SystemExit(
            "No usable cell-type labels were detected; all exported labels would be 'unknown'.\n"
            "Inspect %s and rerun with --cell-type-column <COLUMN>, or set validation.scrna.cell_type_column."
            % (os.path.join(out_dir or ".", "scrna_metadata_columns.tsv"))
        )
    return meta


def downsample_obs(meta: pd.DataFrame, max_cells_per_celltype: Optional[int], seed: int) -> np.ndarray:
    if not max_cells_per_celltype or max_cells_per_celltype <= 0:
        return np.arange(meta.shape[0], dtype=int)
    rng = np.random.default_rng(seed)
    keep = []
    for _, idx in meta.groupby("cell_type", dropna=False).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        if idx.size > max_cells_per_celltype:
            idx = rng.choice(idx, size=max_cells_per_celltype, replace=False)
        keep.extend(idx.tolist())
    return np.asarray(sorted(keep), dtype=int)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out-dir", default="data/external/scrna")
    ap.add_argument("--genes", default="data/external/scrna/scrna_gene_universe.tsv")
    ap.add_argument("--max-cells-per-celltype", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260531)
    ap.add_argument("--cell-type-column", default=None)
    ap.add_argument("--condition-column", default=None)
    ap.add_argument("--sample-column", default=None)
    ap.add_argument("--allow-unknown-celltype", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    genes = read_gene_universe(args.genes)
    if not genes:
        raise SystemExit("No gene universe found. Expected a TSV with a 'gene' column: %s" % args.genes)

    print("Reading h5ad metadata:", args.h5ad)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    symbol_index = make_symbol_index(adata)
    symbol_index = symbol_index[symbol_index["gene"].isin(set(genes))].copy()
    if symbol_index.empty:
        raise SystemExit("None of the requested atlas genes were found in the h5ad var names/metadata.")

    obs = adata.obs.copy().reset_index(drop=True)
    cell_ids = [str(x) for x in list(adata.obs_names)]
    meta_full = infer_metadata(
        obs,
        cell_ids,
        cell_type_column=args.cell_type_column,
        condition_column=args.condition_column,
        sample_column=args.sample_column,
        fail_if_unknown_celltype=not args.allow_unknown_celltype,
        out_dir=args.out_dir,
    )
    print("Detected cell-type labels:", meta_full["cell_type"].nunique(), "unique labels")
    print(meta_full["cell_type"].value_counts(dropna=False).head(20).to_string())
    keep_obs_pos = downsample_obs(meta_full, args.max_cells_per_celltype, args.seed)
    meta = meta_full.iloc[keep_obs_pos, :].copy().reset_index(drop=True)

    var_pos = symbol_index["var_pos"].to_numpy(dtype=int)
    symbols = symbol_index["gene"].tolist()

    print(
        "Exporting compact scRNA table:",
        len(symbols), "genes x", len(keep_obs_pos), "cells",
    )
    X = adata[keep_obs_pos, var_pos].X
    if sp is not None and sp.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)
    X = X.astype(np.float32, copy=False)

    expr = pd.DataFrame(X.T, index=symbols, columns=meta["cell"].astype(str).tolist())
    expr.insert(0, "gene", expr.index)
    expr_path = os.path.join(args.out_dir, "scrna_expression.tsv.gz")
    meta_path = os.path.join(args.out_dir, "scrna_metadata.tsv")
    expr.to_csv(expr_path, sep="\t", index=False, compression="gzip")
    meta.to_csv(meta_path, sep="\t", index=False)

    print("Wrote", expr_path)
    print("Wrote", meta_path)
    print("Python:", sys.executable)


if __name__ == "__main__":
    main()
