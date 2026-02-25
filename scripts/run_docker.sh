#!/usr/bin/env bash
# =============================================================================
# run_docker.sh - Run sc_tools pipeline in Docker for a project
# =============================================================================
# The Docker image uses conda env sc_tools (PATH set in image; no explicit activate).
#
# Usage (from repo root):
#   ./scripts/run_docker.sh <project_path> [command]
#
# Examples:
#   ./scripts/run_docker.sh projects/visium/ggo_visium
#   ./scripts/run_docker.sh projects/visium/ggo_visium python scripts/loupe2adata.py
#   ./scripts/run_docker.sh projects/imc/lymph_dlbcl bash
#
# Or from a project dir:
#   ./run_docker.sh              # if project has run_docker.sh wrapper
#
# Docker is auto-installed on macOS (brew) if missing. On Linux, install manually.
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

IMAGE="${SC_TOOLS_IMAGE:-sc_tools:latest}"

# -----------------------------------------------------------------------------
# Ensure Docker is installed and running
# -----------------------------------------------------------------------------
ensure_docker() {
  if command -v docker &>/dev/null; then
    if docker info &>/dev/null 2>&1; then
      return 0
    fi
    echo "Docker is installed but not running. Start Docker Desktop and try again." >&2
    exit 1
  fi

  echo "Docker not found. Attempting to install..."
  case "$(uname -s)" in
    Darwin)
      if command -v brew &>/dev/null; then
        echo "Running: brew install --cask docker"
        brew install --cask docker
        echo ""
        echo "Docker Desktop installed. Please open Docker Desktop, wait for it to start,"
        echo "then run this script again."
        exit 1
      else
        echo "Install Homebrew (https://brew.sh) first, or install Docker from https://docker.com" >&2
        exit 1
      fi
      ;;
    Linux)
      echo "Install Docker manually:" >&2
      echo "  Ubuntu/Debian: sudo apt-get install docker.io" >&2
      echo "  Fedora:        sudo dnf install docker" >&2
      echo "  Or see: https://docs.docker.com/engine/install/" >&2
      exit 1
      ;;
    *)
      echo "Install Docker from https://docker.com" >&2
      exit 1
      ;;
  esac
}

ensure_docker

# -----------------------------------------------------------------------------
# Build image if missing
# -----------------------------------------------------------------------------
if ! docker image inspect "$IMAGE" &>/dev/null; then
  echo "Image $IMAGE not found. Building..."
  docker build -t "$IMAGE" "$REPO_ROOT"
fi

# -----------------------------------------------------------------------------
# Run container: mount full repo, workdir = project
# Set GGO_VISIUM_PROJECT_DIR when running ggo_visium so score_gene_signatures.py
# resolves paths correctly regardless of which script path is invoked.
# -----------------------------------------------------------------------------
DOCKER_OPTS=(-v "$REPO_ROOT:/workspace" -w "/workspace/$PROJECT_REL")
if [[ "$PROJECT_REL" == "projects/visium/ggo_visium" ]]; then
  DOCKER_OPTS+=(-e "GGO_VISIUM_PROJECT_DIR=/workspace/$PROJECT_REL")
fi
if [[ ${#USER_CMD[@]} -eq 0 ]]; then
  exec docker run -it --rm "${DOCKER_OPTS[@]}" "$IMAGE" bash
else
  # No -t when running a command (e.g. from make) so "not a TTY" does not fail
  exec docker run --rm "${DOCKER_OPTS[@]}" "$IMAGE" "${USER_CMD[@]}"
fi
