#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_ROOT="$ROOT/results/external_deconvolution_inputs/cibersortx"
REF="$INPUT_ROOT/PDAC_7class_refsample.txt"
MIX_DIR="$INPUT_ROOT/mixtures"
RUN_ROOT="$INPUT_ROOT/docker_runs"
RESULT_ROOT="$INPUT_ROOT/results"

: "${CIBERSORTX_EMAIL:?Set CIBERSORTX_EMAIL to your registered academic/non-commercial CIBERSORTx account email.}"
: "${CIBERSORTX_TOKEN:?Set CIBERSORTX_TOKEN to the token supplied by CIBERSORTx.}"

command -v docker >/dev/null 2>&1 || {
  echo "Docker is not available on PATH." >&2
  exit 2
}

[[ -s "$REF" ]] || {
  echo "Missing reference file: $REF" >&2
  echo "Run: Rscript scripts/export_cibersortx_final_validation_inputs.R" >&2
  exit 2
}
[[ -d "$MIX_DIR" ]] || {
  echo "Missing mixture directory: $MIX_DIR" >&2
  exit 2
}

mkdir -p "$RUN_ROOT" "$RESULT_ROOT"

datasets=(GSE15471 GSE28735 GSE62165 GSE16515 GSE71989 GSE91035)

for ds in "${datasets[@]}"; do
  mix="$MIX_DIR/${ds}_mixture.txt"
  [[ -s "$mix" ]] || {
    echo "Missing mixture: $mix" >&2
    exit 2
  }

  in_dir="$RUN_ROOT/${ds}/input"
  out_dir="$RUN_ROOT/${ds}/output"
  rm -rf "$in_dir" "$out_dir"
  mkdir -p "$in_dir" "$out_dir"

  cp "$REF" "$in_dir/refsample.txt"
  cp "$mix" "$in_dir/mixture.txt"

  echo
  echo "========================================================================"
  echo "CIBERSORTx S-mode fractions: $ds"
  echo "========================================================================"

  docker run --rm \
    -v "$in_dir:/src/data:z" \
    -v "$out_dir:/src/outdir:z" \
    cibersortx/fractions \
      --username "$CIBERSORTX_EMAIL" \
      --token "$CIBERSORTX_TOKEN" \
      --single_cell TRUE \
      --refsample refsample.txt \
      --mixture mixture.txt \
      --rmbatchSmode TRUE \
      --QN TRUE \
      --perm 50 \
      --absolute FALSE

  result="$(find "$out_dir" -maxdepth 1 -type f \
    \( -iname '*Adjusted*.txt' -o -iname '*Results*.txt' \) \
    | sort | head -n 1 || true)"

  if [[ -z "$result" || ! -s "$result" ]]; then
    echo "CIBERSORTx did not produce a non-empty fractions result for $ds." >&2
    echo "Output directory: $out_dir" >&2
    find "$out_dir" -maxdepth 1 -type f -printf '  %f\n' >&2 || true
    exit 3
  fi

  final="$RESULT_ROOT/${ds}_CIBERSORTx_Adjusted.txt"
  cp "$result" "$final"
  echo "Saved: $final"
done

echo
echo "All six CIBERSORTx runs completed."
echo "Importing and benchmarking returned fractions..."
cd "$ROOT"
Rscript scripts/import_cibersortx_final_validation.R

echo
echo "Final CIBERSORTx validation complete."
echo "Review:"
echo "  results/tables/external/reference_deconvolution_expected_pair_meta.tsv"
echo "  results/tables/external/reference_deconvolution_qc.tsv"
