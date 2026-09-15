#!/usr/bin/env python3
"""
Create a small TCGA-PAAD / GTEx pancreas expression subset from UCSC Xena Toil.

This is intentionally a stand-alone, line-streaming script, not a targets target.
It avoids the OS-level kills that can happen when R tries to read the full
TcgaTargetGtex_rsem_gene_tpm.gz matrix inside the pipeline.

Inputs:
  - Xena Toil expression matrix: genes x samples, gzipped TSV
  - Xena phenotype table: samples x annotations, gzipped TSV
  - atlas requested gene list: one HUGO symbol per line
  - optional GENCODE v23 probemap for Ensembl -> HUGO symbol mapping

Outputs:
  - data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz
  - data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

TCGA_RE = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", re.I)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace", newline="")
    return open(path, mode, encoding="utf-8", errors="replace", newline="")


def norm_gene(x: str) -> str:
    x = (x or "").strip()
    if "|" in x:
        parts = [p for p in x.split("|") if p]
        # Prefer a non-Ensembl symbol-like part when present.
        for p in reversed(parts):
            if not p.upper().startswith("ENSG"):
                x = p
                break
        else:
            x = parts[-1]
    # Remove Ensembl version suffix only when the identifier looks Ensembl-like.
    if re.match(r"^ENSG\d+\.\d+$", x, re.I):
        x = x.split(".", 1)[0]
    return x.strip().upper()


def strip_ens_version(x: str) -> str:
    x = (x or "").strip()
    if re.match(r"^ENSG\d+(\.\d+)?$", x, re.I):
        return x.split(".", 1)[0].upper()
    return x.upper()


def download_if_missing(url: str, dest: str | Path) -> None:
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    eprint(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)


def read_gene_list(path: str | Path) -> List[str]:
    genes = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            g = norm_gene(line)
            if g:
                genes.append(g)
    genes = sorted(set(genes))
    if len(genes) < 3:
        raise SystemExit(f"Gene list too small: {path}")
    return genes


def read_probemap(path: str | Path | None, wanted_symbols: set[str]) -> Tuple[set[str], Dict[str, str]]:
    """Return wanted row identifiers and id->symbol mapping."""
    wanted_ids: set[str] = set(wanted_symbols)
    id_to_symbol: Dict[str, str] = {}
    if path is None or not Path(path).exists():
        eprint("No probemap found; matching expression row IDs directly to symbols only.")
        return wanted_ids, id_to_symbol

    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return wanted_ids, id_to_symbol
        h = [x.strip().lower() for x in header]
        # Common Xena probemap columns: id, gene, chrom, chromStart, chromEnd, strand
        id_idx = h.index("id") if "id" in h else 0
        if "gene" in h:
            gene_idx = h.index("gene")
        elif "symbol" in h:
            gene_idx = h.index("symbol")
        elif "gene_symbol" in h:
            gene_idx = h.index("gene_symbol")
        else:
            gene_idx = 1 if len(header) > 1 else 0
        n = 0
        matched = 0
        for row in reader:
            if len(row) <= max(id_idx, gene_idx):
                continue
            rid_raw = row[id_idx].strip()
            sym = norm_gene(row[gene_idx])
            rid = strip_ens_version(rid_raw)
            if rid and sym:
                id_to_symbol[rid] = sym
                id_to_symbol[rid_raw.upper()] = sym
            n += 1
            if sym in wanted_symbols:
                wanted_ids.add(rid)
                wanted_ids.add(rid_raw.upper())
                matched += 1
        eprint(f"Loaded probemap rows={n}; atlas symbols with probemap IDs={matched}.")
    return wanted_ids, id_to_symbol


def infer_groups(pheno_path: str | Path, max_samples: int = 1000) -> List[dict]:
    with open_text(pheno_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Empty phenotype file: {pheno_path}")
        names = reader.fieldnames
        lower = {n.lower(): n for n in names}
        sample_col = None
        for cand in ["sample", "sampleid", "sample_id", "_sample", "x.sample", "id"]:
            if cand in lower:
                sample_col = lower[cand]
                break
        if sample_col is None:
            sample_col = names[0]
        out = []
        seen = set()
        for row in reader:
            sample = (row.get(sample_col) or "").strip()
            if not sample or sample in seen:
                continue
            collapsed = " | ".join((v or "") for v in row.values()).lower()
            is_tcga = "tcga" in collapsed or sample.upper().startswith("TCGA")
            is_gtex = "gtex" in collapsed or sample.upper().startswith("GTEX")
            is_pancreas = "pancreas" in collapsed or "pancreatic" in collapsed
            is_paad = (
                "paad" in collapsed
                or "pancreatic adenocarcinoma" in collapsed
                or "pancreatic ductal adenocarcinoma" in collapsed
                or ("pancreas" in collapsed and is_tcga)
            )
            is_primary = (
                "primary tumor" in collapsed
                or "primary solid tumor" in collapsed
                or "tumor" in collapsed
                or bool(TCGA_RE.search(sample))
            )
            condition = cohort = None
            if is_tcga and is_paad and is_primary:
                condition, cohort = "tumor", "TCGA_PAAD"
            elif is_gtex and is_pancreas:
                condition, cohort = "control", "GTEx_pancreas"
            if condition:
                seen.add(sample)
                out.append({"sample": sample, "condition": condition, "cohort": cohort})
    if not out:
        raise SystemExit(
            "No TCGA-PAAD/GTEx pancreas samples were inferred. "
            "Inspect the Xena phenotype file and provide preprocessed metadata if needed."
        )
    if len(out) > max_samples:
        raise SystemExit(
            f"Selected {len(out)} samples, which exceeds max_samples={max_samples}. "
            "The phenotype parser probably selected all TCGA tumors."
        )
    counts = {}
    for r in out:
        k = (r["cohort"], r["condition"])
        counts[k] = counts.get(k, 0) + 1
    eprint("Selected samples: " + "; ".join(f"{a}/{b}={n}" for (a, b), n in sorted(counts.items())))
    return out


def read_expression_header(expr_path: str | Path) -> List[str]:
    with open_text(expr_path) as fh:
        line = fh.readline()
    if not line:
        raise SystemExit(f"Expression file appears empty: {expr_path}")
    return line.rstrip("\n\r").split("\t")


def row_candidates(raw_id: str) -> List[str]:
    raw = (raw_id or "").strip()
    out = {raw.upper(), strip_ens_version(raw), norm_gene(raw)}
    if "|" in raw:
        for p in raw.split("|"):
            out.add(p.strip().upper())
            out.add(strip_ens_version(p))
            out.add(norm_gene(p))
    return [x for x in out if x]


def map_row_to_symbol(raw_id: str, id_to_symbol: Dict[str, str]) -> str:
    for cand in row_candidates(raw_id):
        if cand in id_to_symbol:
            return id_to_symbol[cand]
    return norm_gene(raw_id)


def stream_subset(expr_path: str | Path, metadata: List[dict], wanted_ids: set[str], id_to_symbol: Dict[str, str], out_expr: str | Path) -> None:
    header = read_expression_header(expr_path)
    gene_col = header[0]
    sample_order = [r["sample"] for r in metadata]
    col_index = {name: i for i, name in enumerate(header)}
    found_samples = [s for s in sample_order if s in col_index]
    if len(found_samples) < 10:
        raise SystemExit(f"Too few selected samples found in expression header: {len(found_samples)}")
    if len(found_samples) < len(sample_order):
        missing = len(sample_order) - len(found_samples)
        eprint(f"Warning: {missing} selected phenotype samples were not in the expression header.")
    idx = [0] + [col_index[s] for s in found_samples]
    max_idx = max(idx)

    Path(out_expr).parent.mkdir(parents=True, exist_ok=True)
    total = kept = 0
    with open_text(expr_path) as fh, gzip.open(out_expr, "wt", encoding="utf-8", newline="") as out:
        _ = fh.readline()
        out.write("gene\t" + "\t".join(found_samples) + "\n")
        for line in fh:
            total += 1
            if total % 5000 == 0:
                eprint(f"Scanned {total} expression rows; retained {kept}...", end="\r")
            raw = line.split("\t", 1)[0].strip()
            if not any(c in wanted_ids for c in row_candidates(raw)):
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) <= max_idx:
                continue
            symbol = map_row_to_symbol(raw, id_to_symbol)
            values = [parts[i] for i in idx[1:]]
            out.write(symbol + "\t" + "\t".join(values) + "\n")
            kept += 1
    eprint(f"\nFinished expression subset: scanned rows={total}; retained rows={kept}.")
    if kept < 3:
        raise SystemExit(
            "Retained fewer than 3 expression rows. The gene identifiers probably do not match. "
            "Ensure the GENCODE v23 probemap is present or provide a preprocessed expression file."
        )


def write_metadata(metadata: List[dict], out_meta: str | Path, expression_header_samples: Sequence[str] | None = None) -> None:
    Path(out_meta).parent.mkdir(parents=True, exist_ok=True)
    if expression_header_samples is not None:
        keep = set(expression_header_samples)
        metadata = [r for r in metadata if r["sample"] in keep]
    with open(out_meta, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "condition", "cohort"], delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--expression", default="data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz")
    ap.add_argument("--phenotype", default="data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz")
    ap.add_argument("--gene-list", default="data/external/cache/tcga_gtex_requested_genes.txt")
    ap.add_argument("--probemap", default="data/external/tcga_gtex/gencode.v23.annotation.gene.probemap")
    ap.add_argument("--probemap-url", default="https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap")
    ap.add_argument("--out-expression", default="data/external/tcga_gtex/tcga_paad_gtex_pancreas_expression.tsv.gz")
    ap.add_argument("--out-metadata", default="data/external/tcga_gtex/tcga_paad_gtex_pancreas_metadata.tsv")
    ap.add_argument("--max-samples", type=int, default=1000)
    ap.add_argument("--no-download-probemap", action="store_true")
    args = ap.parse_args(argv)

    for required in [args.expression, args.phenotype, args.gene_list]:
        if not Path(required).exists():
            raise SystemExit(f"Required file not found: {required}")

    if not args.no_download_probemap and args.probemap_url:
        download_if_missing(args.probemap_url, args.probemap)

    genes = read_gene_list(args.gene_list)
    eprint(f"Requested atlas genes: {len(genes)}")
    wanted_ids, id_to_symbol = read_probemap(args.probemap, set(genes))
    metadata = infer_groups(args.phenotype, max_samples=args.max_samples)
    stream_subset(args.expression, metadata, wanted_ids, id_to_symbol, args.out_expression)
    # Keep metadata in inferred order. R will filter to samples present in expression if needed.
    write_metadata(metadata, args.out_metadata)
    eprint(f"Wrote expression: {args.out_expression}")
    eprint(f"Wrote metadata:   {args.out_metadata}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
