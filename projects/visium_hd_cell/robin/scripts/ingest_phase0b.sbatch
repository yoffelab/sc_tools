#!/bin/bash
#SBATCH --job-name=robin_ingest_cellseg
#SBATCH --partition=scu-cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=/athena/elementolab/scratch/juk4007/robin/sc_tools/projects/visium_hd_cell/robin/logs/ingest_cellseg_%j.log

# Phase 0b: Ingest Visium HD cell segmentation data for robin project
# Needs more memory due to geojson parsing (large polygon files)

set -euo pipefail

cd /athena/elementolab/scratch/juk4007/robin/sc_tools

# Activate conda
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

# Ensure sc_tools is importable
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"

echo "=== Phase 0b: Visium HD cell segmentation ingest ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

python projects/visium_hd_cell/robin/scripts/ingest_phase0b.py

echo "End: $(date)"
