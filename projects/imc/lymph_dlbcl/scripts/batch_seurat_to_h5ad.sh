#!/usr/bin/env bash
# Batch convert Seurat RDS files under a directory to AnnData h5ad.
# Run from project root: projects/imc/lymph_dlbcl/
#
# Usage:
#   ./scripts/batch_seurat_to_h5ad.sh [rds_root]
#   ./scripts/batch_seurat_to_h5ad.sh [rds_root] --list
#
# Arguments:
#   rds_root  Directory containing .rds files (default: data/downloaded).
#   --list    If set, only convert files listed in metadata/files_to_download_rds.csv
#             (relative_path column); paths are relative to rds_root.
#
# Output: results/seurat_converted/<relative_path>.h5ad for each .rds

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

USE_LIST=""
RDS_ROOT="data/downloaded"
if [[ "${1:-}" == "--list" ]]; then
  USE_LIST=1
elif [[ "${2:-}" == "--list" ]]; then
  USE_LIST=1
  RDS_ROOT="${1}"
else
  [[ -n "${1:-}" ]] && RDS_ROOT="$1"
fi

if [[ ! -d "$RDS_ROOT" ]]; then
  echo "Error: RDS root directory not found: $RDS_ROOT" >&2
  exit 1
fi

R_SCRIPT="$SCRIPT_DIR/seurat_rds_to_h5ad.R"
if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: R script not found: $R_SCRIPT" >&2
  exit 1
fi

OUT_ROOT="results/seurat_converted"
mkdir -p "$OUT_ROOT"

if [[ -n "$USE_LIST" ]]; then
  LIST_FILE="metadata/files_to_download_rds.csv"
  if [[ ! -f "$LIST_FILE" ]]; then
    echo "Error: List file not found: $LIST_FILE" >&2
    exit 1
  fi
  # Skip header; use first column (relative_path)
  tail -n +2 "$LIST_FILE" | cut -d',' -f1 | while IFS= read -r rel; do
    rel="${rel%"${rel##*[![:space:]]}"}"
    [[ -z "$rel" ]] && continue
    in_path="$RDS_ROOT/$rel"
    if [[ ! -f "$in_path" ]]; then
      echo "Skip (not found): $in_path" >&2
      continue
    fi
    out_rel="${rel%.rds}.h5ad"
    out_path="$OUT_ROOT/$out_rel"
    mkdir -p "$(dirname "$out_path")"
    echo "Converting: $rel -> $out_rel"
    Rscript "$R_SCRIPT" "$in_path" "$out_path" "$PROJECT_ROOT" || true
  done
else
  find "$RDS_ROOT" -type f -name '*.rds' | while IFS= read -r in_path; do
    rel="${in_path#$RDS_ROOT/}"
    out_rel="${rel%.rds}.h5ad"
    out_path="$OUT_ROOT/$out_rel"
    mkdir -p "$(dirname "$out_path")"
    echo "Converting: $rel -> $out_rel"
    Rscript "$R_SCRIPT" "$in_path" "$out_path" "$PROJECT_ROOT" || true
  done
fi

echo "Done. Outputs under $OUT_ROOT/"
