#!/bin/bash
#SBATCH --job-name=scanvi_retry
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --array=0,2,31,32,33,34,35,36,37,38,39
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/scanvi_rsearch_retry_%A_%a.log
#
# Retry failed/timed-out tasks from job 2703212.
# Tasks 0,2: failed due to /tmp race condition (fixed: task-specific temp path)
# Tasks 31-39: timed out at 2h (fixed: 4h limit)
#
# Submit:
#   sbatch scripts/run_scanvi_rsearch_retry.sh

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial

mkdir -p ${WORKDIR}/results/scanvi_rsearch
mkdir -p ${WORKDIR}/logs

echo "============================================"
echo "scANVI retry: task ${SLURM_ARRAY_TASK_ID} of 40"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Start: $(date)"
echo "============================================"

# Fix: use task-specific temp path to avoid race condition
INPUT=${WORKDIR}/results/m2_benchmark/adata.m2.h5ad
LOCAL_INPUT=/tmp/adata_rsearch_${SLURM_ARRAY_TASK_ID}.h5ad

if [ ! -f "${LOCAL_INPUT}" ]; then
    echo "Copying input to ${LOCAL_INPUT} for faster I/O..."
    cp ${INPUT} ${LOCAL_INPUT}
    echo "Done copying ($(du -h ${LOCAL_INPUT} | cut -f1))"
else
    echo "Input already cached at ${LOCAL_INPUT} ($(du -h ${LOCAL_INPUT} | cut -f1))"
fi

# Override INPUT_PATH env var so the python script uses the local copy
export SCANVI_INPUT_OVERRIDE=${LOCAL_INPUT}

# Run the worker
python ${WORKDIR}/scripts/run_scanvi_rsearch.py

# Cleanup temp file
rm -f ${LOCAL_INPUT}

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
