#!/bin/bash
#SBATCH --job-name=m2_benchmark
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/m2_final_benchmark_%j.log
#
# M2 final benchmark: train top 3 scANVI configs, compute UMAPs for all
# methods, generate comparison report.
#
# Estimated: ~60-90 min (3 scANVI trainings + 10 UMAP computations + plotting)
#
# Submit:
#   sbatch scripts/run_m2_final_benchmark.sh

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial
cd ${WORKDIR}

echo "============================================"
echo "M2 Final Benchmark"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Start: $(date)"
echo "============================================"

python scripts/run_m2_final_benchmark.py

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
