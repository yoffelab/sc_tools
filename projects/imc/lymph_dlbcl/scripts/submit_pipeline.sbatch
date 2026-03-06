#!/bin/bash
#SBATCH --job-name=dlbcl_imc_pipeline
#SBATCH --partition=scu-cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=outputs/logs/pipeline_%j.out
#SBATCH --error=outputs/logs/pipeline_%j.err

set -eo pipefail

echo "=== DLBCL IMC Pipeline ==="
echo "Job ID: ${SLURM_JOB_ID:-nohup}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo ""

cd /home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl

# Activate environment
source ~/.bashrc
conda activate sc_tools 2>/dev/null || true

echo "Python: $(which python)"
echo "Working dir: $(pwd)"
echo ""

# Ensure output dirs exist
mkdir -p outputs/logs figures/manuscript figures/QC metadata results/seurat_converted

# h5ad files already uploaded to results/seurat_converted/
echo "=== Checking h5ad files ==="
ls -lh results/seurat_converted/tcell_2_preprocessing/t2_SO_seurat.h5ad 2>/dev/null || echo "WARNING: immune h5ad missing"
ls -lh results/seurat_converted/stroma_2_preprocessing/S1_seurat_SO.h5ad 2>/dev/null || echo "WARNING: stromal h5ad missing"
echo ""

# ============================================================
# Phase 0: Download clinical metadata + validate h5ad
# ============================================================
echo "=== Phase 0: Download + Validate ==="
bash scripts/download_clinical_metadata.sh 2>&1 || echo "WARNING: Some downloads may have failed (check manifest)"
python scripts/validate_h5ad_objects.py 2>&1

# ============================================================
# Phase 1: Build AnnData checkpoints
# ============================================================
echo ""
echo "=== Phase 1: Build AnnData ==="
python scripts/build_panel_adata.py --panel both 2>&1
python scripts/attach_clinical_metadata.py --panel both 2>&1
python scripts/build_spatial_adata.py 2>&1

# ============================================================
# Phase 2-3: Validate cell types + assign LME
# ============================================================
echo ""
echo "=== Phase 2-3: Validate + LME ==="
python scripts/validate_celltypes.py --panel both 2>&1
python scripts/build_lme_classes.py --panel both 2>&1

# ============================================================
# Phase 4: Main figures
# ============================================================
echo ""
echo "=== Phase 4: Main Figures ==="
python scripts/fig1_single_cell_atlas.py 2>&1 || echo "WARNING: Fig 1 failed"
python scripts/fig2_lme_classes.py 2>&1 || echo "WARNING: Fig 2 failed"
python scripts/fig3_clinical.py 2>&1 || echo "WARNING: Fig 3 failed"
python scripts/fig4_spatial.py 2>&1 || echo "WARNING: Fig 4 failed"
python scripts/fig5_ml_framework.py 2>&1 || echo "WARNING: Fig 5 failed"

# Supplementary figures
echo ""
echo "=== Phase 4: Supplementary Figures ==="
python scripts/supp_fig1_qc_panels.py 2>&1 || echo "WARNING: Supp Fig 1 failed"
python scripts/supp_fig2_bcell.py 2>&1 || echo "WARNING: Supp Fig 2 failed"
python scripts/supp_fig3_tcell_myeloid.py 2>&1 || echo "WARNING: Supp Fig 3 failed"
python scripts/supp_fig4_vessel.py 2>&1 || echo "WARNING: Supp Fig 4 failed"
python scripts/supp_fig5_tme_sensitivity.py 2>&1 || echo "WARNING: Supp Fig 5 failed"
python scripts/supp_fig6_mutations.py 2>&1 || echo "WARNING: Supp Fig 6 failed"
python scripts/supp_fig7_rna_protein.py 2>&1 || echo "WARNING: Supp Fig 7 failed"
python scripts/supp_fig8_extended_survival.py 2>&1 || echo "WARNING: Supp Fig 8 failed"

# ============================================================
# Phase 6: Validation
# ============================================================
echo ""
echo "=== Phase 6: Validation ==="
python scripts/validate_figures.py 2>&1

echo ""
echo "=== Pipeline complete ==="
echo "Finished: $(date)"
echo "Check outputs/ and figures/ for results"

# List generated figures
echo ""
echo "=== Generated figures ==="
find figures/ -name "*.pdf" -type f 2>/dev/null | sort
