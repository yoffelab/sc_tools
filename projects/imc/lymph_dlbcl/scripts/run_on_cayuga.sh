#!/usr/bin/env bash
# run_on_cayuga.sh — Sync project to cayuga and run Snakemake pipeline
#
# Usage:
#   bash scripts/run_on_cayuga.sh [phase0|phase1|phase2|phase3|figures|all]
#   bash scripts/run_on_cayuga.sh phase0   # Just download + validate
#   bash scripts/run_on_cayuga.sh all      # Full pipeline
#
# Prerequisites:
#   - SSH access to cayuga configured (ssh cayuga works)
#   - Data present at /home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl/

set -euo pipefail

REMOTE_HOST="cayuga"
REMOTE_PROJECT="/home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl"
LOCAL_PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
TARGET="${1:-all}"

echo "=== DLBCL IMC Pipeline — cayuga ==="
echo "Local:  ${LOCAL_PROJECT}"
echo "Remote: ${REMOTE_HOST}:${REMOTE_PROJECT}"
echo "Target: ${TARGET}"
echo ""

# Step 1: Sync scripts, config, and Snakefile to cayuga
echo "--- Syncing project files to cayuga ---"
rsync -avz --exclude='*.h5ad' --exclude='data/' --exclude='notebooks/' \
    --exclude='results/seurat_converted/' --exclude='manuscript/' \
    --exclude='.DS_Store' --exclude='__pycache__/' \
    "${LOCAL_PROJECT}/" "${REMOTE_HOST}:${REMOTE_PROJECT}/"

echo ""
echo "--- Running Snakemake on cayuga ---"

# Step 2: Run on cayuga
# Use conda activate or module load as needed
ssh "${REMOTE_HOST}" bash -l << REMOTE_SCRIPT
set -euo pipefail
cd "${REMOTE_PROJECT}"

# Activate environment (adjust as needed for your HPC setup)
if command -v conda &>/dev/null; then
    conda activate sc_tools 2>/dev/null || true
fi

echo "Working directory: \$(pwd)"
echo "Python: \$(which python)"
echo ""

# Ensure output dirs exist
mkdir -p outputs figures/manuscript figures/QC metadata

# Run Snakemake
if command -v snakemake &>/dev/null; then
    snakemake --cores 8 ${TARGET} --printshellcmds 2>&1 | tee outputs/snakemake_${TARGET}_\$(date +%Y%m%d_%H%M%S).log
else
    echo "snakemake not found; running scripts directly..."

    case "${TARGET}" in
        phase0)
            echo "=== Phase 0: Download + Validate ==="
            bash scripts/download_clinical_metadata.sh
            python scripts/validate_h5ad_objects.py
            ;;
        phase1)
            echo "=== Phase 1: Build AnnData ==="
            python scripts/build_panel_adata.py --panel both
            python scripts/attach_clinical_metadata.py --panel both
            python scripts/build_spatial_adata.py
            ;;
        phase3)
            echo "=== Phase 3: Validate + LME ==="
            python scripts/validate_celltypes.py --panel both
            python scripts/build_lme_classes.py --panel both
            ;;
        figures|fig1|fig2|fig3|fig4|fig5)
            echo "=== Phase 4: Figures ==="
            python scripts/fig1_single_cell_atlas.py
            python scripts/fig2_lme_classes.py
            python scripts/fig3_clinical.py
            python scripts/fig4_spatial.py
            python scripts/fig5_ml_framework.py
            for n in 1 2 3 4 5 6 7 8; do
                python scripts/supp_fig\${n}_*.py 2>&1 || true
            done
            ;;
        all)
            echo "=== Full Pipeline ==="
            bash scripts/download_clinical_metadata.sh
            python scripts/validate_h5ad_objects.py
            python scripts/build_panel_adata.py --panel both
            python scripts/attach_clinical_metadata.py --panel both
            python scripts/build_spatial_adata.py
            python scripts/validate_celltypes.py --panel both
            python scripts/build_lme_classes.py --panel both
            python scripts/fig1_single_cell_atlas.py
            python scripts/fig2_lme_classes.py
            python scripts/fig3_clinical.py
            python scripts/fig4_spatial.py
            python scripts/fig5_ml_framework.py
            for n in 1 2 3 4 5 6 7 8; do
                python scripts/supp_fig\${n}_*.py 2>&1 || true
            done
            python scripts/validate_figures.py
            ;;
        *)
            echo "Unknown target: ${TARGET}"
            echo "Usage: run_on_cayuga.sh [phase0|phase1|phase3|figures|all]"
            exit 1
            ;;
    esac
fi

echo ""
echo "=== Done ==="
REMOTE_SCRIPT

echo ""
echo "--- Syncing results back ---"
# Pull back outputs and figures
rsync -avz "${REMOTE_HOST}:${REMOTE_PROJECT}/outputs/" "${LOCAL_PROJECT}/outputs/"
rsync -avz "${REMOTE_HOST}:${REMOTE_PROJECT}/figures/" "${LOCAL_PROJECT}/figures/"
rsync -avz "${REMOTE_HOST}:${REMOTE_PROJECT}/metadata/lme_class_assignments.csv" "${LOCAL_PROJECT}/metadata/" 2>/dev/null || true
rsync -avz "${REMOTE_HOST}:${REMOTE_PROJECT}/metadata/download_manifest_phase0.csv" "${LOCAL_PROJECT}/metadata/" 2>/dev/null || true

echo ""
echo "=== Pipeline complete ==="
echo "Check outputs/ and figures/ for results"
