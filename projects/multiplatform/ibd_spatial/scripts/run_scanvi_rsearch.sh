#!/bin/bash
#SBATCH --job-name=scanvi_rsearch
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --array=0-39%8
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/scanvi_rsearch_%A_%a.log
#
# scANVI random hyperparameter search — 40 configs on cayuga (A100 80GB)
#
# Array: 0-39, max 8 concurrent (polite to other users)
# Each task: ~20-40 min on A100 (164K cells, 2000 HVGs)
# Total wall time: ~40 tasks / 8 concurrent x 40 min = ~3.3h
#
# Submit:
#   sbatch scripts/run_scanvi_rsearch.sh
#
# Monitor:
#   squeue -u juk4007 -j <jobid>
#   sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,NodeList
#
# After completion:
#   python scripts/aggregate_scanvi_rsearch.py

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial

# Create output and log directories
mkdir -p ${WORKDIR}/results/scanvi_rsearch
mkdir -p ${WORKDIR}/logs

echo "============================================"
echo "scANVI random search: task ${SLURM_ARRAY_TASK_ID} of 40"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Start: $(date)"
echo "============================================"

# Copy input data to /tmp for faster I/O on Lustre
INPUT=${WORKDIR}/results/m2_benchmark/adata.m2.h5ad
LOCAL_INPUT=/tmp/adata.m2.h5ad

if [ ! -f "${LOCAL_INPUT}" ]; then
    echo "Copying input to /tmp for faster I/O..."
    cp ${INPUT} ${LOCAL_INPUT}
    echo "Done copying ($(du -h ${LOCAL_INPUT} | cut -f1))"
else
    echo "Input already cached at ${LOCAL_INPUT} ($(du -h ${LOCAL_INPUT} | cut -f1))"
fi

# Run the worker (reads SLURM_ARRAY_TASK_ID from env)
python ${WORKDIR}/scripts/run_scanvi_rsearch.py

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
