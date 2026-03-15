#!/bin/bash
#SBATCH --job-name=resolvi_ss
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --nodelist=g0001
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/resolvi_semisup.log

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial

python ${WORKDIR}/scripts/run_resolvi_semisup.py
