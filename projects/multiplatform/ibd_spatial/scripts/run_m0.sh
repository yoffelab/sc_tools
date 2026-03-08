#!/bin/bash
#SBATCH --job-name=m0_bench
#SBATCH --partition=scu-cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/m0_benchmark.log

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial

python ${WORKDIR}/scripts/run_m0_benchmark.py
