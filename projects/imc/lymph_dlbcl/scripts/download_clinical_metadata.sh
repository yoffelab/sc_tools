#!/usr/bin/env bash
# download_clinical_metadata.sh
# Phase 0.1: Download clinical/metadata files from cayuga backup
#
# Usage: bash scripts/download_clinical_metadata.sh
# Run from: projects/imc/lymph_dlbcl/
#
# Supports running ON cayuga (cp) or locally (scp from cayuga)

set -eo pipefail

BACKUP_DLBCLV2="/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"
BACKUP_META="/home/fs01/juk4007/elementolab/backup/dylan/data/hyperion/DLBCL"
LOCAL_CLINICAL="data/downloaded/clinical"
LOCAL_METADATA="data/downloaded/metadata"
MANIFEST="metadata/download_manifest_phase0.csv"

mkdir -p "$LOCAL_CLINICAL" "$LOCAL_METADATA"

# Initialize manifest
echo "remote_path,local_path,category,description,download_date" > "$MANIFEST"

# Detect if running on cayuga (local copy) vs remote (scp)
if [ -d "$BACKUP_DLBCLV2" ]; then
    ON_CAYUGA=true
    echo "Running on cayuga — using local cp"
else
    ON_CAYUGA=false
    echo "Running remotely — using scp from cayuga"
fi

download_file() {
    local remote_path="$1"
    local local_dir="$2"
    local category="$3"
    local description="$4"
    local filename
    filename=$(basename "$remote_path")
    local local_path="${local_dir}/${filename}"

    echo "  ${filename}"
    if [ "$ON_CAYUGA" = true ]; then
        if cp "$remote_path" "$local_path" 2>/dev/null; then
            echo "${remote_path},${local_path},${category},${description},$(date +%Y-%m-%d)" >> "$MANIFEST"
            return 0
        fi
    else
        if scp "cayuga:${remote_path}" "$local_path" 2>/dev/null; then
            echo "${remote_path},${local_path},${category},${description},$(date +%Y-%m-%d)" >> "$MANIFEST"
            return 0
        fi
    fi
    echo "    -> NOT FOUND"
    return 1
}

echo "=== Phase 0.1: Downloading clinical and metadata files ==="
echo ""

# --- DLBCLv2 root files ---
echo "--- DLBCLv2 root ---"
download_file "$BACKUP_DLBCLV2/2.17.22.v2stroma_C8_clusters.csv" "$LOCAL_METADATA" "tme" "TME C8 cluster assignments" || true
download_file "$BACKUP_DLBCLV2/2.23.23.TME_alt6.csv" "$LOCAL_METADATA" "tme" "TME alt 6-class assignments" || true
download_file "$BACKUP_DLBCLV2/2.23.23.IHC_337_cases.csv" "$LOCAL_METADATA" "ihc" "IHC data 337 cases" || true
download_file "$BACKUP_DLBCLV2/2.23.23.IHC_266_all_markers.csv" "$LOCAL_METADATA" "ihc" "IHC all markers 266 cases" || true
download_file "$BACKUP_DLBCLV2/s1_tme_table.csv" "$LOCAL_METADATA" "tme" "Slide 1 TME table" || true
download_file "$BACKUP_DLBCLV2/s2_tme_table.csv" "$LOCAL_METADATA" "tme" "Slide 2 TME table" || true
download_file "$BACKUP_DLBCLV2/s_tme_clusters.csv" "$LOCAL_METADATA" "tme" "TME clusters" || true
download_file "$BACKUP_DLBCLV2/s_tme_all.csv" "$LOCAL_METADATA" "tme" "TME all data" || true
download_file "$BACKUP_DLBCLV2/s1_cellID_table.csv" "$LOCAL_METADATA" "cellid" "S1 cell ID table" || true
download_file "$BACKUP_DLBCLV2/s2_cellID_table.csv" "$LOCAL_METADATA" "cellid" "S2 cell ID table" || true
download_file "$BACKUP_DLBCLV2/full_coo_c.csv" "$LOCAL_METADATA" "clinical" "Full COO classification" || true
download_file "$BACKUP_DLBCLV2/percent_meta.csv" "$LOCAL_METADATA" "clinical" "Percent metadata" || true
download_file "$BACKUP_DLBCLV2/staudt_clinical_cyto_clusters.csv" "$LOCAL_METADATA" "clinical" "Staudt clinical cyto clusters" || true

# --- clinical_analysis/ ---
echo ""
echo "--- clinical_analysis ---"
download_file "$BACKUP_DLBCLV2/clinical_analysis/7.14.22.TME.zscore.csv" "$LOCAL_METADATA" "tme" "TME z-scores" || true
download_file "$BACKUP_DLBCLV2/clinical_analysis/7.14.22.totals.risk.csv" "$LOCAL_METADATA" "clinical" "Survival risk classification" || true
download_file "$BACKUP_DLBCLV2/clinical_analysis/8.15.22.DLC380_patient_tme.csv" "$LOCAL_METADATA" "tme" "Patient-level TME (380 cases)" || true
download_file "$BACKUP_DLBCLV2/clinical_analysis/4.14.22.merged_covariates_clin.csv" "$LOCAL_CLINICAL" "clinical" "Merged clinical covariates" || true
download_file "$BACKUP_DLBCLV2/clinical_analysis/2.25.23.TME_prediction_test.csv" "$LOCAL_METADATA" "ml" "TME prediction test data" || true
download_file "$BACKUP_DLBCLV2/clinical_analysis/7.18.22.stringent_pos_abundance.csv" "$LOCAL_METADATA" "abundance" "Stringent positive abundance" || true

# --- total_tumor/ ---
echo ""
echo "--- total_tumor ---"
download_file "$BACKUP_DLBCLV2/total_tumor/1.13.21.merged.abundance.csv" "$LOCAL_METADATA" "abundance" "Merged abundance table" || true
download_file "$BACKUP_DLBCLV2/total_tumor/merged_abundance.csv" "$LOCAL_METADATA" "abundance" "Merged abundance (alt)" || true
download_file "$BACKUP_DLBCLV2/total_tumor/5.27.22.merged.counts.csv" "$LOCAL_METADATA" "abundance" "Merged cell counts" || true

# --- stroma_spatial/ ---
echo ""
echo "--- stroma_spatial ---"
download_file "$BACKUP_DLBCLV2/stroma_spatial/6.21.22.community_merged_sc_df_stroma.csv" "$LOCAL_METADATA" "spatial" "Community merged single-cell stroma" || true
download_file "$BACKUP_DLBCLV2/stroma_spatial/7.31.22.k30_per_patient_community_abundance.csv" "$LOCAL_METADATA" "spatial" "Per-patient community abundance k30" || true
download_file "$BACKUP_DLBCLV2/stroma_spatial/3.18.23.KNN30_communities_abun.csv" "$LOCAL_METADATA" "spatial" "KNN30 community abundance" || true

# --- Cell ID CSVs ---
echo ""
echo "--- Cell ID mapping ---"
download_file "$BACKUP_DLBCLV2/stroma_1_preprocessing/S1_seurat_cellid.csv" "$LOCAL_METADATA" "cellid" "S1 Seurat cell ID" || true
download_file "$BACKUP_DLBCLV2/stroma_2_preprocessing/S2_seurat_cellid.csv" "$LOCAL_METADATA" "cellid" "S2 Seurat cell ID" || true
download_file "$BACKUP_DLBCLV2/tcell_2_preprocessing/t2_seurat_cellid.csv" "$LOCAL_METADATA" "cellid" "T2 Seurat cell ID" || true

# --- mutations/ ---
echo ""
echo "--- mutations ---"
download_file "$BACKUP_DLBCLV2/mutations/8.14.22.CK.deconv.table.csv" "$LOCAL_METADATA" "mutations" "Mutation deconv table" || true

# --- Clinical from DLBCL meta (different backup path) ---
echo ""
echo "--- DLBCL meta (alternate paths) ---"
download_file "$BACKUP_META/meta/BCCA/DLBCL_clinical_full.csv" "$LOCAL_CLINICAL" "clinical" "Full clinical data" || true
download_file "$BACKUP_META/meta_data/CTMA_121_punch_notes.csv" "$LOCAL_CLINICAL" "clinical" "Sample-to-slide mapping" || true
download_file "$BACKUP_META/meta_data/CTMA121_mut_table.csv" "$LOCAL_CLINICAL" "clinical" "Mutation table per gene" || true

echo ""
echo "=== Download complete ==="
echo "Manifest: ${MANIFEST}"
NFILES=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
echo "Downloaded: ${NFILES} files"
echo ""
if [ "$NFILES" -lt 10 ]; then
    echo "WARNING: Fewer than 10 files downloaded. Some paths may need fixing."
fi
