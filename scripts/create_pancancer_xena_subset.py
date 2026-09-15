#!/usr/bin/env python3
"""
Create a memory-safe pan-cancer specificity subset from the UCSC Xena Toil
TCGA/TARGET/GTEx expression matrix already used by the PDAC validation.

Panel:
  TCGA COAD tumors vs GTEx colon
  TCGA STAD tumors vs GTEx stomach
  TCGA LIHC tumors vs GTEx liver
  TCGA LUAD tumors vs GTEx lung

Only the current atlas signature-universe genes are retained.

Outputs:
  data/external/pancancer/pancancer_expression.tsv.gz
  data/external/pancancer/pancancer_metadata.tsv
  data/external/pancancer/pancancer_cohort_counts.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

TCGA_PRIMARY_RE = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", re.I)

PANEL = {
    "COAD": {
        "disease_tokens": ("COAD", "COLON ADENOCARCINOMA", "COLON CANCER"),
        "gtex_site": "COLON",
    },
    "STAD": {
        "disease_tokens": ("STAD", "STOMACH ADENOCARCINOMA", "GASTRIC ADENOCARCINOMA"),
        "gtex_site": "STOMACH",
    },
    "LIHC": {
        "disease_tokens": ("LIHC", "LIVER HEPATOCELLULAR CARCINOMA", "HEPATOCELLULAR CARCINOMA"),
        "gtex_site": "LIVER",
    },
    "LUAD": {
        "disease_tokens": ("LUAD", "LUNG ADENOCARCINOMA"),
        "gtex_site": "LUNG",
    },
}


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
        for p in reversed(parts):
            if not p.upper().startswith("ENSG"):
                x = p
                break
        else:
            x = parts[-1]
    if re.match(r"^ENSG\d+\.\d+$", x, re.I):
        x = x.split(".", 1)[0]
    return x.strip().upper()


def strip_ens_version(x: str) -> str:
    x = (x or "").strip()
    if re.match(r"^ENSG\d+(\.\d+)?$", x, re.I):
        return x.split(".", 1)[0].upper()
    return x.upper()


def read_gene_list(path: str | Path) -> List[str]:
    out = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            g = norm_gene(line)
            if g:
                out.append(g)
    out = sorted(set(out))
    if len(out) < 3:
        raise SystemExit(f"Gene list too small: {path}")
    return out


def read_probemap(path: str | Path, wanted_symbols: set[str]) -> Tuple[set[str], Dict[str, str]]:
    wanted_ids = set(wanted_symbols)
    id_to_symbol: Dict[str, str] = {}
    if not Path(path).exists():
        eprint("Probemap absent; direct symbol matching only.")
        return wanted_ids, id_to_symbol

    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        h = [x.strip().lower() for x in header]
        id_idx = h.index("id") if "id" in h else 0
        if "gene" in h:
            gene_idx = h.index("gene")
        elif "symbol" in h:
            gene_idx = h.index("symbol")
        elif "gene_symbol" in h:
            gene_idx = h.index("gene_symbol")
        else:
            gene_idx = 1 if len(header) > 1 else 0
        matched = n = 0
        for row in reader:
            if len(row) <= max(id_idx, gene_idx):
                continue
            rid_raw = row[id_idx].strip()
            rid = strip_ens_version(rid_raw)
            sym = norm_gene(row[gene_idx])
            if rid and sym:
                id_to_symbol[rid] = sym
                id_to_symbol[rid_raw.upper()] = sym
            n += 1
            if sym in wanted_symbols:
                wanted_ids.add(rid)
                wanted_ids.add(rid_raw.upper())
                matched += 1
        eprint(f"Loaded probemap rows={n}; requested symbols with mapped IDs={matched}.")
    return wanted_ids, id_to_symbol


def _norm_colname(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (x or "").strip().lower()).strip()


def first_existing_col(fieldnames: List[str], candidates: Sequence[str]) -> str | None:
    # Exact case-insensitive match first, then punctuation/whitespace-normalized
    # matching. UCSC Xena phenotype releases have used both underscore-style
    # and human-readable headers (e.g. "primary disease or tissue").
    lower = {x.lower(): x for x in fieldnames}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    normalized = {_norm_colname(x): x for x in fieldnames}
    for c in candidates:
        nc = _norm_colname(c)
        if nc in normalized:
            return normalized[nc]
    return None


def site_matches(site: str, wanted: str) -> bool:
    s = (site or "").upper()
    w = wanted.upper()
    if w == "COLON":
        return "COLON" in s
    return w in s


def infer_groups(pheno_path: str | Path) -> List[dict]:
    """
    Infer the pre-specified TCGA/GTEx specificity panel from the combined
    UCSC Xena Toil phenotype file.

    Important: in TcgaTargetGTEX_phenotype.txt, `_study` identifies the broad
    source (TCGA/TARGET/GTEX), while the cancer diagnosis is held in the
    disease / `_primary_disease` field. Do not use `_study` to infer COAD,
    STAD, LIHC or LUAD.
    """
    with open_text(pheno_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Empty phenotype file: {pheno_path}")

        f = reader.fieldnames
        sample_col = first_existing_col(
            f, ["sample", "sampleid", "sample_id", "_sample", "id"]
        ) or f[0]
        study_col = first_existing_col(f, ["_study", "study"])
        disease_col = first_existing_col(
            f,
            [
                "_primary_disease",
                "primary_disease",
                "primary disease",
                "primary disease or tissue",
                "disease",
                "cancer_type",
                "cancer type abbreviation",
                "project_id",
                "project",
            ],
        )
        site_col = first_existing_col(f, ["_primary_site", "primary_site"])
        sample_type_col = first_existing_col(f, ["_sample_type", "sample_type"])
        detail_col = first_existing_col(
            f, ["body_site_detail (SMTSD)", "body_site_detail", "SMTSD"]
        )
        category_col = first_existing_col(
            f, ["_category", "category", "detailed_category", "TCGA_GTEX_main_category"]
        )

        if study_col is None:
            raise SystemExit(
                "Could not find _study/study in Xena phenotype table. "
                f"Columns were: {', '.join(f)}"
            )
        if disease_col is None:
            # Last-resort semantic header detection for future Xena releases.
            disease_like = [
                col for col in f
                if "disease" in _norm_colname(col)
                or ("cancer" in _norm_colname(col) and "type" in _norm_colname(col))
            ]
            if len(disease_like) == 1:
                disease_col = disease_like[0]
            else:
                raise SystemExit(
                    "Could not find a unique disease/cancer column in Xena phenotype table. "
                    f"Columns were: {', '.join(f)}"
                )

        eprint(
            "Phenotype mapping: "
            f"sample={sample_col}; study={study_col}; disease={disease_col}; "
            f"primary_site={site_col}; sample_type={sample_type_col}; category={category_col}"
        )

        out = []
        seen = set()
        tcga_primary_n = 0
        tcga_disease_counts: Dict[str, int] = {}

        for row in reader:
            sample = (row.get(sample_col) or "").strip()
            if not sample or sample in seen:
                continue

            study = (row.get(study_col) or "").strip()
            disease = (row.get(disease_col) or "").strip()
            site = (row.get(site_col) or "").strip() if site_col else ""
            detail = (row.get(detail_col) or "").strip() if detail_col else ""
            category = (row.get(category_col) or "").strip() if category_col else ""
            stype = (row.get(sample_type_col) or "").strip() if sample_type_col else ""

            study_u = study.upper()
            disease_u = disease.upper()
            site_u = (site or detail).upper()
            category_u = category.upper()
            stype_u = stype.upper()
            sample_u = sample.upper()

            is_gtex = "GTEX" in study_u or sample_u.startswith("GTEX")
            is_tcga = "TCGA" in study_u or sample_u.startswith("TCGA")

            if is_tcga:
                is_primary = (
                    ("PRIMARY" in stype_u and "TUMOR" in stype_u)
                    or ("PRIMARY SOLID TUMOR" in stype_u)
                    or bool(TCGA_PRIMARY_RE.search(sample_u))
                )
                if not is_primary:
                    continue

                tcga_primary_n += 1
                if disease:
                    tcga_disease_counts[disease] = tcga_disease_counts.get(disease, 0) + 1

                # Primary classification: disease / primary disease.
                # Secondary fallbacks use category/study only when they contain
                # an explicit cancer token; primary site alone is deliberately
                # NOT used for lung, because it would mix LUAD and LUSC.
                combined = " | ".join([disease_u, category_u, study_u])
                matched = None
                for cancer, spec in PANEL.items():
                    if any(tok in combined for tok in spec["disease_tokens"]):
                        matched = cancer
                        break

                if matched is not None:
                    out.append(
                        {
                            "sample": sample,
                            "cancer_type": matched,
                            "class": "tumor",
                            "study": study,
                            "disease": disease,
                            "primary_site": site or detail,
                            "sample_type": stype,
                        }
                    )
                    seen.add(sample)

            elif is_gtex:
                effective_site = site if site else detail
                for cancer, spec in PANEL.items():
                    if site_matches(effective_site, spec["gtex_site"]):
                        out.append(
                            {
                                "sample": sample,
                                "cancer_type": cancer,
                                "class": "gtex_normal",
                                "study": study,
                                "disease": disease,
                                "primary_site": effective_site,
                                "sample_type": stype,
                            }
                        )
                        seen.add(sample)
                        break

    counts = {}
    for r in out:
        k = (r["cancer_type"], r["class"])
        counts[k] = counts.get(k, 0) + 1

    eprint(
        "Selected samples: "
        + "; ".join(f"{c}/{cl}={n}" for (c, cl), n in sorted(counts.items()))
    )

    missing_groups = []
    for cancer in PANEL:
        for cl in ("tumor", "gtex_normal"):
            if counts.get((cancer, cl), 0) < 20:
                missing_groups.append(
                    f"{cancer}/{cl}={counts.get((cancer, cl), 0)}"
                )

    if missing_groups:
        top_diseases = sorted(
            tcga_disease_counts.items(), key=lambda kv: (-kv[1], kv[0])
        )[:40]
        diag = "; ".join(f"{k}={v}" for k, v in top_diseases)
        raise SystemExit(
            "Pan-cancer panel has fewer than 20 samples in required groups: "
            + ", ".join(missing_groups)
            + f". TCGA primary-tumor rows seen={tcga_primary_n}. "
            + "Top TCGA disease labels: "
            + (diag if diag else "<none>")
        )

    return out

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


def read_expression_header(expr_path: str | Path) -> List[str]:
    with open_text(expr_path) as fh:
        line = fh.readline()
    if not line:
        raise SystemExit(f"Expression file appears empty: {expr_path}")
    return line.rstrip("\n\r").split("\t")


def stream_subset(expr_path: str | Path, metadata: List[dict], wanted_ids: set[str],
                  id_to_symbol: Dict[str, str], out_expr: str | Path) -> List[str]:
    header = read_expression_header(expr_path)
    sample_order = [r["sample"] for r in metadata]
    col_index = {name: i for i, name in enumerate(header)}
    found = [s for s in sample_order if s in col_index]
    if len(found) < 100:
        raise SystemExit(f"Too few selected samples in expression matrix: {len(found)}")
    if len(found) < len(sample_order):
        eprint(f"Warning: {len(sample_order) - len(found)} phenotype samples lack expression columns.")

    idx = [col_index[s] for s in found]
    max_idx = max(idx)
    Path(out_expr).parent.mkdir(parents=True, exist_ok=True)

    total = kept = 0
    with open_text(expr_path) as fh, gzip.open(out_expr, "wt", encoding="utf-8", newline="") as out:
        _ = fh.readline()
        out.write("gene\t" + "\t".join(found) + "\n")
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
            out.write(symbol + "\t" + "\t".join(parts[i] for i in idx) + "\n")
            kept += 1

    eprint(f"\nFinished pan-cancer subset: scanned rows={total}; retained rows={kept}; samples={len(found)}.")
    if kept < 100:
        raise SystemExit("Too few signature-universe genes retained; check probemap/gene identifiers.")
    return found


def write_metadata(metadata: List[dict], found: Sequence[str], out_meta: str | Path, out_counts: str | Path) -> None:
    keep = set(found)
    rows = [r for r in metadata if r["sample"] in keep]
    Path(out_meta).parent.mkdir(parents=True, exist_ok=True)
    with open(out_meta, "w", encoding="utf-8", newline="") as fh:
        fields = ["sample", "cancer_type", "class", "study", "disease", "primary_site", "sample_type"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    counts = {}
    for r in rows:
        k = (r["cancer_type"], r["class"])
        counts[k] = counts.get(k, 0) + 1
    with open(out_counts, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cancer_type", "class", "n"])
        for (c, cl), n in sorted(counts.items()):
            w.writerow([c, cl, n])


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--expression", default="data/external/tcga_gtex/TcgaTargetGtex_rsem_gene_tpm.gz")
    ap.add_argument("--phenotype", default="data/external/tcga_gtex/TcgaTargetGTEX_phenotype.txt.gz")
    ap.add_argument("--gene-list", default="data/external/cache/tcga_gtex_requested_genes.txt")
    ap.add_argument("--probemap", default="data/external/tcga_gtex/gencode.v23.annotation.gene.probemap")
    ap.add_argument("--out-expression", default="data/external/pancancer/pancancer_expression.tsv.gz")
    ap.add_argument("--out-metadata", default="data/external/pancancer/pancancer_metadata.tsv")
    ap.add_argument("--out-counts", default="data/external/pancancer/pancancer_cohort_counts.tsv")
    args = ap.parse_args(argv)

    for p in [args.expression, args.phenotype, args.gene_list]:
        if not Path(p).exists():
            raise SystemExit(f"Required local Xena input not found: {p}")

    genes = read_gene_list(args.gene_list)
    eprint(f"Requested atlas genes: {len(genes)}")
    wanted_ids, id_to_symbol = read_probemap(args.probemap, set(genes))
    metadata = infer_groups(args.phenotype)
    found = stream_subset(args.expression, metadata, wanted_ids, id_to_symbol, args.out_expression)
    write_metadata(metadata, found, args.out_metadata, args.out_counts)
    eprint(f"Wrote expression: {args.out_expression}")
    eprint(f"Wrote metadata:   {args.out_metadata}")
    eprint(f"Wrote counts:     {args.out_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
