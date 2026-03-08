#!/bin/bash
#SBATCH --job-name=robin_cell_pipeline
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=/athena/elementolab/scratch/juk4007/robin/sc_tools/projects/visium_hd_cell/robin/logs/pipeline_%j.log

# Robin cell segmentation pipeline: qc_filter -> metadata_attach -> preprocess -> scoring
# 1.57M cells x 18K genes — needs GPU for scVI in integration benchmark
# Preprocess runs full integration benchmark: Harmony, BBKNN, ComBat, Scanorama, scVI, PCA

set -euo pipefail

cd /athena/elementolab/scratch/juk4007/robin/sc_tools

# Activate conda
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

# Ensure sc_tools is importable
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"

PROJECT="projects/visium_hd_cell/robin"

echo "=== Robin Cell Pipeline ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "Memory: $(free -h | head -2)"
nvidia-smi 2>/dev/null || echo "No GPU detected (CPU fallback for scVI)"

# Phase 1: QC filter
echo ""
echo "=== Phase 1: QC Filter ==="
echo "Start: $(date)"
python ${PROJECT}/scripts/qc_filter.py
echo "Phase 1 done: $(date)"

# Phase 2: Metadata attachment
echo ""
echo "=== Phase 2: Metadata Attach ==="
echo "Start: $(date)"
python ${PROJECT}/scripts/attach_metadata.py
echo "Phase 2 done: $(date)"

# Phase 3: Preprocessing with integration benchmark
echo ""
echo "=== Phase 3: Preprocess (Integration Benchmark) ==="
echo "Start: $(date)"
python ${PROJECT}/scripts/preprocess.py
echo "Phase 3 done: $(date)"

# Phase 3.5b: Gene signature scoring
echo ""
echo "=== Phase 3.5b: Scoring ==="
echo "Start: $(date)"
python ${PROJECT}/scripts/score_signatures.py
echo "Phase 3.5b done: $(date)"

echo ""
echo "=== Pipeline complete ==="
echo "End: $(date)"
