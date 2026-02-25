#!/usr/bin/env bash
# =============================================================================
# run_apptainer.sh - Run sc_tools pipeline in Apptainer for a project
# =============================================================================
# Uses Apptainer (Singularity) for portability (HPC, no Docker daemon).
# SIF is built from the repo Dockerfile; see project_setup.md.
#
# Usage (from repo root):
#   ./scripts/run_apptainer.sh <project_path> [command]
#
# Examples:
#   ./scripts/run_apptainer.sh projects/visium/ggo_visium
#   ./scripts/run_apptainer.sh projects/visium/ggo_visium python scripts/loupe2adata.py
#   ./scripts/run_apptainer.sh projects/visium_hd/robin python scripts/run_signature_scoring.py
#
# Requires: Apptainer (or Singularity) installed; SIF at containers/sc_tools.sif
# Build SIF: see project_setup.md (e.g. apptainer build from docker-daemon).
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_PATH="${1:?Usage: $0 <project_path> [command...]}"
shift || true
USER_CMD=("$@")

# Resolve project path (relative to repo root)
if [[ "$PROJECT_PATH" != /* ]]; then
  PROJECT_PATH="$REPO_ROOT/$PROJECT_PATH"
fi
PROJECT_PATH="$(cd "$PROJECT_PATH" && pwd)"
PROJECT_REL="${PROJECT_PATH#$REPO_ROOT/}"

if [[ ! "$PROJECT_PATH" == "$REPO_ROOT"/* ]]; then
  echo "Error: Project must be under repo root: $REPO_ROOT" >&2
  exit 1
fi

SIF="${SC_TOOLS_SIF:-$REPO_ROOT/containers/sc_tools.sif}"

# -----------------------------------------------------------------------------
# Ensure Apptainer (or Singularity) and SIF exist
# -----------------------------------------------------------------------------
if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
  echo "Error: Apptainer or Singularity not found. Install Apptainer (apptainer.org) or Singularity." >&2
  exit 1
fi

APPTAINER_CMD=""
if command -v apptainer &>/dev/null; then
  APPTAINER_CMD=apptainer
elif command -v singularity &>/dev/null; then
  APPTAINER_CMD=singularity
fi

if [[ ! -f "$SIF" ]]; then
  echo "Error: SIF not found: $SIF" >&2
  echo "Build it with: apptainer build $SIF docker-daemon://sc_tools:latest" >&2
  echo "(after: docker build -t sc_tools:latest .)" >&2
  echo "See project_setup.md for details." >&2
  exit 1
fi

# -----------------------------------------------------------------------------
# Run: bind repo to /workspace, workdir = project; no -it when command given
# -----------------------------------------------------------------------------
BIND="$REPO_ROOT:/workspace"
WORKDIR="/workspace/$PROJECT_REL"

# ggo_visium: set env so scripts resolve project dir
if [[ "$PROJECT_REL" == "projects/visium/ggo_visium" ]]; then
  export APPTAINERENV_GGO_VISIUM_PROJECT_DIR="$WORKDIR"
fi

if [[ ${#USER_CMD[@]} -eq 0 ]]; then
  exec "$APPTAINER_CMD" exec --bind "$BIND" --pwd "$WORKDIR" "$SIF" bash
else
  exec "$APPTAINER_CMD" exec --bind "$BIND" --pwd "$WORKDIR" "$SIF" "${USER_CMD[@]}"
fi
