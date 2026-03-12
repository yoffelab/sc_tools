#!/bin/bash
#SBATCH --job-name=resolvi_scanvi
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a40:1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/resolvi_scanvi.log

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial

python ${WORKDIR}/scripts/run_resolvi_scanvi.py
