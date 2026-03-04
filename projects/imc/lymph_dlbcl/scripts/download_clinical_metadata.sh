#!/usr/bin/env bash
# download_clinical_metadata.sh
# Phase 0.1: Download clinical/metadata files from cayuga backup
#
# Usage: bash scripts/download_clinical_metadata.sh
# Run from: projects/imc/lymph_dlbcl/
#
# Remote base: /home/fs01/juk4007/elementolab/backup/dylan/
# Local dest:  data/downloaded/clinical/ and data/downloaded/metadata/

set -euo pipefail

REMOTE_HOST="cayuga"
REMOTE_BASE="/home/fs01/juk4007/elementolab/backup/dylan"
LOCAL_CLINICAL="data/downloaded/clinical"
LOCAL_METADATA="data/downloaded/metadata"
MANIFEST="metadata/download_manifest_phase0.csv"

mkdir -p "$LOCAL_CLINICAL" "$LOCAL_METADATA"

# Initialize manifest
echo "remote_path,local_path,category,description,download_date" > "$MANIFEST"

download_file() {
    local remote_path="$1"
    local local_dir="$2"
    local category="$3"
    local description="$4"
    local filename
    filename=$(basename "$remote_path")
    local local_path="${local_dir}/${filename}"

    echo "Downloading: ${remote_path}"
    if scp "${REMOTE_HOST}:${REMOTE_BASE}/${remote_path}" "$local_path" 2>/dev/null; then
        echo "${REMOTE_BASE}/${remote_path},${local_path},${category},${description},$(date +%Y-%m-%d)" >> "$MANIFEST"
        echo "  -> OK: ${local_path}"
    else
        echo "  -> FAILED: ${remote_path} (may not exist at this path)"
        # Try alternative paths
        return 1
    fi
}

echo "=== Phase 0.1: Downloading clinical and metadata files ==="
echo ""

# --- Clinical data (Fig 3, 5) ---
echo "--- Clinical data ---"

download_file \
    "data/hyperion/DLBCL/meta/BCCA/DLBCL_clinical_full.csv" \
    "$LOCAL_CLINICAL" "clinical" "Full clinical data (OS/PFS/COO/LymphGen)" || \
download_file \
    "hyperion/DLBCLv2/total_tumor/DLBCL_clinical_full.csv" \
    "$LOCAL_CLINICAL" "clinical" "Full clinical data (alt path)" || \
    echo "  WARNING: DLBCL_clinical_full.csv not found at known paths"

download_file \
    "data/hyperion/DLBCL/meta_data/CTMA_121_punch_notes.csv" \
    "$LOCAL_CLINICAL" "clinical" "Sample-to-slide mapping" || \
    echo "  WARNING: CTMA_121_punch_notes.csv not found"

download_file \
    "data/hyperion/DLBCL/meta_data/CTMA121_mut_table.csv" \
    "$LOCAL_CLINICAL" "clinical" "Mutation table per gene" || \
    echo "  WARNING: CTMA121_mut_table.csv not found"

# --- TME class assignments (Fig 2) ---
echo ""
echo "--- TME class assignments ---"

download_file \
    "hyperion/DLBCLv2/2.17.22.v2stroma_C8_clusters.csv" \
    "$LOCAL_METADATA" "tme" "TME/stromal C8 cluster assignments" || \
    echo "  WARNING: v2stroma_C8_clusters.csv not found"

download_file \
    "hyperion/DLBCLv2/2.23.23.TME_alt6.csv" \
    "$LOCAL_METADATA" "tme" "Alternative TME 6-class assignments" || \
    echo "  WARNING: TME_alt6.csv not found"

download_file \
    "hyperion/DLBCLv2/7.14.22.TME.zscore.csv" \
    "$LOCAL_METADATA" "tme" "TME z-scores per sample" || \
    echo "  WARNING: TME.zscore.csv not found"

# --- IHC data (Fig 5) ---
echo ""
echo "--- IHC data ---"

download_file \
    "hyperion/DLBCLv2/2.23.23.IHC_337_cases.csv" \
    "$LOCAL_METADATA" "ihc" "IHC data for 337 cases" || \
    echo "  WARNING: IHC_337_cases.csv not found"

download_file \
    "hyperion/DLBCLv2/2.23.23.IHC_266_all_markers.csv" \
    "$LOCAL_METADATA" "ihc" "IHC all markers for 266 cases" || \
    echo "  WARNING: IHC_266_all_markers.csv not found"

# --- Cell abundance tables (Fig 1, 2) ---
echo ""
echo "--- Cell abundance ---"

download_file \
    "hyperion/DLBCLv2/total_tumor/1.13.21.merged.abundance.csv" \
    "$LOCAL_METADATA" "abundance" "Merged cell type abundance table" || \
    echo "  WARNING: merged.abundance.csv not found"

# --- Spatial community data (Fig 4) ---
echo ""
echo "--- Spatial community ---"

download_file \
    "hyperion/DLBCLv2/stroma_spatial/6.21.22.community_merged_sc_df_stroma.csv" \
    "$LOCAL_METADATA" "spatial" "Spatial community merged single-cell data" || \
    echo "  WARNING: community_merged_sc_df_stroma.csv not found"

# --- Cell ID CSVs (for sample mapping) ---
echo ""
echo "--- Cell ID mapping ---"

download_file \
    "hyperion/DLBCLv2/stroma_2_preprocessing/S2_seurat_cellid.csv" \
    "$LOCAL_METADATA" "cellid" "Stromal S2 cell ID mapping" || \
    echo "  WARNING: S2_seurat_cellid.csv not found"

download_file \
    "hyperion/DLBCLv2/stroma_1_preprocessing/S1_seurat_cellid.csv" \
    "$LOCAL_METADATA" "cellid" "Stromal S1 cell ID mapping" || \
    echo "  WARNING: S1_seurat_cellid.csv not found"

# --- TME tables ---
echo ""
echo "--- TME tables ---"

download_file \
    "hyperion/DLBCLv2/total_tumor/s1_tme_table.csv" \
    "$LOCAL_METADATA" "tme" "Slide 1 TME abundance table" || \
    echo "  WARNING: s1_tme_table.csv not found"

download_file \
    "hyperion/DLBCLv2/total_tumor/s2_tme_table.csv" \
    "$LOCAL_METADATA" "tme" "Slide 2 TME abundance table" || \
    echo "  WARNING: s2_tme_table.csv not found"

# --- Patient-level TME ---
download_file \
    "hyperion/DLBCLv2/8.15.22.DLC380_patient_tme.csv" \
    "$LOCAL_METADATA" "tme" "Patient-level TME table (380 cases)" || \
    echo "  WARNING: DLC380_patient_tme.csv not found"

# --- Survival/risk ---
download_file \
    "hyperion/DLBCLv2/total_tumor/7.14.22.totals.risk.csv" \
    "$LOCAL_METADATA" "clinical" "Survival risk classification table" || \
    echo "  WARNING: totals.risk.csv not found"

# --- Community abundance (Fig 4) ---
download_file \
    "hyperion/DLBCLv2/stroma_spatial/7.31.22.k30_per_patient_community_abundance.csv" \
    "$LOCAL_METADATA" "spatial" "Per-patient community abundance (k=30)" || \
    echo "  WARNING: k30_per_patient_community_abundance.csv not found"

echo ""
echo "=== Download complete ==="
echo "Manifest: ${MANIFEST}"
echo ""
echo "Downloaded files:"
wc -l < "$MANIFEST" | xargs -I{} echo "  {} entries (including header)"
echo ""
echo "Check manifest for any missing files and download manually if needed."
