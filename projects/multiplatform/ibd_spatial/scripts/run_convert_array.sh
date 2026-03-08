#!/bin/bash
#SBATCH --job-name=rds2h5ad
#SBATCH --partition=scu-cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/convert_%a.log
#SBATCH --array=0-51

# Convert all 52 Seurat v5 RDS files to h5ad
# Usage: sbatch run_convert_array.sh
#
# Task mapping (52 total):
#   0-15:  CosMx_1k (16 samples)
#   16-19: CosMx_6k (4 samples)
#   20-35: Xenium_MT (16 samples)
#   36-39: Xenium_5K (4 samples)
#   40-43: Xenium_MT_withseg (4 samples)
#   44-47: Xenium_MT_noseg (4 samples)
#   48-51: Xenium_colon (4 samples)

set -eo pipefail

# Load R module (Rscript not in PATH on compute nodes by default)
# Note: cannot use `set -u` because /etc/profile.d/lang.sh uses unset LC_ALL
source /etc/profile || true
module load R/4.4.1
export LD_LIBRARY_PATH=/opt/ohpc/pub/software/R/4.4.1/lib64:${LD_LIBRARY_PATH:-}
RSCRIPT=/opt/ohpc/pub/software/R/4.4.1/bin/Rscript

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial
SCRIPT=${WORKDIR}/scripts/convert_rds_to_h5ad.R
DATA_IBD=/athena/project-saha/data_IBD
CSV=${DATA_IBD}/NC_compare_11122025_meta.csv
TASK=${SLURM_ARRAY_TASK_ID}

# Ensure R libraries are available (libs_R441 built against R/4.4.1 module)
# Use ONLY libs_R441 to avoid loading .so files compiled against system R
export R_LIBS_USER=~/R/libs_R441

# Ensure Python (for h5ad assembly step) is available
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate base 2>/dev/null || true

# Map task ID to panel, sample number, RDS path, output path
if [ $TASK -le 15 ]; then
    PANEL=cosmx_1k
    SAMPLE_NUM=$((TASK + 1))
    RDS=${DATA_IBD}/CosMx_1k_16/${SAMPLE_NUM}_smi_so_qcnorm.RDS
elif [ $TASK -le 19 ]; then
    PANEL=cosmx_6k
    SAMPLE_NUM=$((TASK - 16 + 1))
    RDS=${DATA_IBD}/CosMx_6k_4/${SAMPLE_NUM}_smi_so_qcnorm.RDS
elif [ $TASK -le 35 ]; then
    PANEL=xenium_mt
    SAMPLE_NUM=$((TASK - 20 + 1))
    RDS=${DATA_IBD}/Xenium_MT_16/${SAMPLE_NUM}_xen_so_qcnorm.RDS
elif [ $TASK -le 39 ]; then
    PANEL=xenium_5k
    SAMPLE_NUM=$((TASK - 36 + 1))
    RDS=${DATA_IBD}/Xenium_5K_4/${SAMPLE_NUM}_xen_so_qcnorm.RDS
elif [ $TASK -le 43 ]; then
    PANEL=xenium_withseg
    SAMPLE_NUM=$((TASK - 40 + 1))
    RDS=${DATA_IBD}/Xenium_MT_withseg_4/${SAMPLE_NUM}_xen_377_with_seg_so_qcnorm.RDS
elif [ $TASK -le 47 ]; then
    PANEL=xenium_noseg
    SAMPLE_NUM=$((TASK - 44 + 1))
    RDS=${DATA_IBD}/Xenium_MT_noseg_4/${SAMPLE_NUM}_xen_377_no_seg_so_qcnorm.RDS
elif [ $TASK -le 51 ]; then
    PANEL=xenium_colon
    SAMPLE_NUM=$((TASK - 48 + 1))
    RDS=${DATA_IBD}/Xenium_colon_4/${SAMPLE_NUM}_xen_colon_with_seg_so_qcnorm.RDS
else
    echo "ERROR: invalid task ID $TASK"
    exit 1
fi

SAMPLE_ID=$(printf "%s_%02d" "$PANEL" "$SAMPLE_NUM")
OUTPUT=${WORKDIR}/data/${SAMPLE_ID}/adata.p0.h5ad

echo "Task $TASK: panel=$PANEL sample=$SAMPLE_NUM rds=$RDS output=$OUTPUT"

# Skip if already converted
if [ -f "$OUTPUT" ]; then
    echo "Output exists, skipping: $OUTPUT"
    exit 0
fi

$RSCRIPT --vanilla "$SCRIPT" "$RDS" "$OUTPUT" "$SAMPLE_NUM" "$PANEL" "$CSV"

echo "Task $TASK completed: $OUTPUT"
