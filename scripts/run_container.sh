#!/usr/bin/env bash
# =============================================================================
# run_container.sh — Unified container runner for sc_tools
# =============================================================================
# Detects OS and available runtime, then runs in the appropriate container:
#   macOS (Darwin) → Docker  (Apptainer requires Linux kernel; macOS = no native)
#   Linux          → Apptainer/Singularity preferred; Docker fallback if missing
#   none           → run directly in active conda env (no container)
#
# Usage (from repo root):
#   ./scripts/run_container.sh <project_path> [command...]
#
# Examples:
#   ./scripts/run_container.sh projects/visium/ggo_visium
#   ./scripts/run_container.sh projects/visium/ggo_visium python scripts/score_gene_signatures.py
#   ./scripts/run_container.sh projects/imc/ggo-imc bash
#
# Environment overrides:
#   SC_TOOLS_RUNTIME=docker|apptainer|singularity|none
#                       Force a specific runtime. Use 'none' to run directly
#                       in the active conda env (no container). Auto-detected
#                       if unset: macOS→docker, Linux→apptainer>singularity>docker>none.
#   SC_TOOLS_IMAGE=myimage:tag           Docker image override (default: sc_tools:latest)
#   SC_TOOLS_SIF=/path/to/sc_tools.sif  SIF path override (default: containers/sc_tools.sif)
#
# All projects get PROJECT_DIR set to the project workdir (inside container or local).
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
PROJECT_PATH="$(cd "$PROJECT_PATH" 2>/dev/null && pwd)" || {
  echo "Error: Project path does not exist: $1" >&2; exit 1
}
if [[ "$PROJECT_PATH" == "$REPO_ROOT" ]]; then
  PROJECT_REL="."
else
  PROJECT_REL="${PROJECT_PATH#"$REPO_ROOT"/}"
  if [[ "$PROJECT_REL" == "$PROJECT_PATH" ]]; then
    echo "Error: Project must be under repo root: $REPO_ROOT" >&2; exit 1
  fi
fi

# =============================================================================
# Runtime detection
# =============================================================================

OS="$(uname -s)"

detect_runtime() {
  if [[ -n "$SC_TOOLS_RUNTIME" ]]; then
    echo "$SC_TOOLS_RUNTIME"
    return
  fi
  case "$OS" in
    Darwin)
      if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
        echo "docker"
      else
        echo "none"
      fi
      ;;
    Linux)
      if command -v apptainer &>/dev/null; then
        echo "apptainer"
      elif command -v singularity &>/dev/null; then
        echo "singularity"
      elif command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
        echo "docker"
      else
        echo "none"
      fi
      ;;
    *)
      if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
        echo "docker"
      else
        echo "none"
      fi
      ;;
  esac
}

RUNTIME="$(detect_runtime)"

echo "[run_container] OS=$OS  RUNTIME=$RUNTIME  PROJECT=$PROJECT_REL"

# =============================================================================
# Docker runner
# =============================================================================

run_docker() {
  local IMAGE="${SC_TOOLS_IMAGE:-sc_tools:latest}"

  # Ensure Docker is running
  if ! command -v docker &>/dev/null; then
    echo "Error: Docker not found." >&2
    if [[ "$OS" == "Darwin" ]] && command -v brew &>/dev/null; then
      echo "Run: brew install --cask docker" >&2
    fi
    exit 1
  fi
  if ! docker info &>/dev/null 2>&1; then
    echo "Error: Docker daemon not running. Start Docker Desktop and retry." >&2
    exit 1
  fi

  # Auto-build image if missing
  if ! docker image inspect "$IMAGE" &>/dev/null 2>&1; then
    echo "[run_container] Image $IMAGE not found. Building from $REPO_ROOT/Dockerfile..."
    docker build -t "$IMAGE" "$REPO_ROOT"
  fi

  local OPTS=()
  if [[ "$PROJECT_REL" == "." || "$PROJECT_REL" == sc_tools* ]]; then
    # sc_tools development: full repo access
    OPTS=(-v "$REPO_ROOT:/workspace" -w "/workspace"
          -e "PROJECT_DIR=/workspace")
  else
    # Project execution: scoped mounts (project RW, sc_tools + scripts RO)
    OPTS=(
      -v "$PROJECT_PATH:/workspace/project"
      -v "$REPO_ROOT/sc_tools:/workspace/sc_tools:ro"
      -v "$REPO_ROOT/scripts:/workspace/scripts:ro"
      -v "$REPO_ROOT/pyproject.toml:/workspace/pyproject.toml:ro"
      -v "$REPO_ROOT/skills.md:/workspace/skills.md:ro"
      -w "/workspace/project"
      -e "PROJECT_DIR=/workspace/project"
    )
  fi

  if [[ ${#USER_CMD[@]} -eq 0 ]]; then
    exec docker run -it --rm "${OPTS[@]}" "$IMAGE" bash
  else
    exec docker run --rm "${OPTS[@]}" "$IMAGE" "${USER_CMD[@]}"
  fi
}

# =============================================================================
# Apptainer/Singularity runner
# =============================================================================

run_apptainer() {
  local SIF="${SC_TOOLS_SIF:-$REPO_ROOT/containers/sc_tools.sif}"
  local CMD="$RUNTIME"  # apptainer or singularity

  if [[ ! -f "$SIF" ]]; then
    echo "Error: SIF not found: $SIF" >&2
    echo "Build it with:" >&2
    echo "  docker build -t sc_tools:latest $REPO_ROOT" >&2
    echo "  apptainer build $SIF docker-daemon://sc_tools:latest" >&2
    echo "Or override: SC_TOOLS_SIF=/path/to/sc_tools.sif $0 ..." >&2
    exit 1
  fi

  local BIND WORKDIR
  if [[ "$PROJECT_REL" == "." || "$PROJECT_REL" == sc_tools* ]]; then
    # sc_tools development: full repo access
    BIND="$REPO_ROOT:/workspace"
    WORKDIR="/workspace"
  else
    # Project execution: scoped mounts (project RW, sc_tools + scripts RO)
    BIND="$PROJECT_PATH:/workspace/project,$REPO_ROOT/sc_tools:/workspace/sc_tools:ro,$REPO_ROOT/scripts:/workspace/scripts:ro,$REPO_ROOT/pyproject.toml:/workspace/pyproject.toml:ro,$REPO_ROOT/skills.md:/workspace/skills.md:ro"
    WORKDIR="/workspace/project"
  fi

  export APPTAINERENV_PROJECT_DIR="$WORKDIR"
  export SINGULARITYENV_PROJECT_DIR="$WORKDIR"  # singularity compat

  if [[ ${#USER_CMD[@]} -eq 0 ]]; then
    exec "$CMD" exec --bind "$BIND" --pwd "$WORKDIR" "$SIF" bash
  else
    exec "$CMD" exec --bind "$BIND" --pwd "$WORKDIR" "$SIF" "${USER_CMD[@]}"
  fi
}

# =============================================================================
# Direct runner (no container — uses active conda env)
# =============================================================================

run_direct() {
  local ACTIVE_ENV="${CONDA_DEFAULT_ENV:-$(basename "$VIRTUAL_ENV" 2>/dev/null)}"
  echo "[run_container] No container runtime — running directly in env: ${ACTIVE_ENV:-<system python>}"
  export PROJECT_DIR="$PROJECT_PATH"

  if [[ ${#USER_CMD[@]} -eq 0 ]]; then
    exec bash --norc -i
  else
    cd "$PROJECT_PATH"
    exec "${USER_CMD[@]}"
  fi
}

# =============================================================================
# Dispatch
# =============================================================================

case "$RUNTIME" in
  docker)        run_docker ;;
  apptainer)     run_apptainer ;;
  singularity)   run_apptainer ;;
  none)          run_direct ;;
  *)
    echo "Error: Unknown runtime: $RUNTIME (valid: docker|apptainer|singularity|none)" >&2
    exit 1
    ;;
esac
