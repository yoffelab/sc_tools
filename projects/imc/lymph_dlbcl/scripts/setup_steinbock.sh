#!/usr/bin/env bash
# setup_steinbock.sh — Pull steinbock Docker image as Apptainer SIF on cayuga.
#
# Run once on cayuga before starting Plan B IMC reprocessing.
#
# Usage:
#   ssh cayuga
#   cd /home/fs01/juk4007/elementolab/sc_tools/projects/imc/lymph_dlbcl
#   bash scripts/setup_steinbock.sh
#
# The SIF is written to ~/elementolab/scratch/steinbock.sif (matches
# config.yaml: phase0a.steinbock_sif).
# Pulling from Docker Hub requires internet access — cayuga login nodes have
# outbound access; compute nodes may not.
#
# Requirements:
#   - Apptainer v1.1+ (module load apptainer or available as apptainer on PATH)
#   - ~5 GB free disk space in ~/elementolab/scratch/
#   - Network access to ghcr.io (GitHub Container Registry)

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
STEINBOCK_IMAGE="ghcr.io/bodenmillergroup/steinbock:latest-cellpose"
SIF_PATH="${HOME}/elementolab/scratch/steinbock.sif"
SCRATCH_DIR="$(dirname "$SIF_PATH")"

# ── Checks ───────────────────────────────────────────────────────────────────
echo "=== steinbock Apptainer setup ==="
echo "Image : $STEINBOCK_IMAGE"
echo "SIF   : $SIF_PATH"
echo ""

# Ensure scratch directory exists and has space
mkdir -p "$SCRATCH_DIR"
AVAILABLE_KB=$(df -k "$SCRATCH_DIR" | awk 'NR==2 {print $4}')
REQUIRED_KB=$((6 * 1024 * 1024))   # 6 GB conservative estimate
if [ "$AVAILABLE_KB" -lt "$REQUIRED_KB" ]; then
    echo "ERROR: not enough free space in $SCRATCH_DIR"
    echo "  Available: $((AVAILABLE_KB / 1024 / 1024)) GB"
    echo "  Required : ~6 GB"
    exit 1
fi

# Load Apptainer module if needed (cayuga uses Lmod)
if ! command -v apptainer &>/dev/null; then
    if command -v module &>/dev/null; then
        module load apptainer 2>/dev/null || true
    fi
fi

if ! command -v apptainer &>/dev/null; then
    echo "ERROR: apptainer not found. Load with: module load apptainer"
    exit 1
fi

APPTAINER_VERSION=$(apptainer --version 2>&1 | head -1)
echo "Apptainer: $APPTAINER_VERSION"

# ── Pull ─────────────────────────────────────────────────────────────────────
if [ -f "$SIF_PATH" ]; then
    echo "SIF already exists: $SIF_PATH"
    echo "Delete and re-run to refresh: rm $SIF_PATH"
    exit 0
fi

echo ""
echo "Pulling $STEINBOCK_IMAGE ..."
echo "(this may take 5-15 minutes on first pull)"
echo ""

# APPTAINER_TMPDIR controls where layers are cached during build
# Use scratch to avoid filling home quota
export APPTAINER_TMPDIR="${SCRATCH_DIR}/apptainer_tmp"
mkdir -p "$APPTAINER_TMPDIR"

apptainer pull "$SIF_PATH" "docker://${STEINBOCK_IMAGE}"

# Verify the SIF runs
echo ""
echo "Verifying steinbock ..."
apptainer exec "$SIF_PATH" steinbock --version

echo ""
echo "=== Setup complete ==="
echo "steinbock SIF: $SIF_PATH"
echo ""
echo "Next steps:"
echo "  1. Confirm raw data paths on cayuga:"
echo "     ls /athena/elementolab/scratch/dym2001/data/hyperion/DLBCL/"
echo "  2. Update metadata/phase0/batch1_immune.tsv and batch1_stromal.tsv"
echo "     with real mcd_file paths (replace placeholder DLBCLv2 paths)."
echo "  3. Enable Plan B in config.yaml: phase0a.enabled: true"
echo "  4. Dry-run: snakemake -n phase0a_immune --config 'phase0a={enabled: true}'"
echo "  5. Submit: snakemake phase0a_immune --executor slurm --jobs 20"
