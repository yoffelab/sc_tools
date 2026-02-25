#!/usr/bin/env bash
# Run Docker for ggo_visium. From repo root: ./projects/visium/ggo_visium/run_docker.sh [command]
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
exec "$REPO_ROOT/scripts/run_docker.sh" projects/visium/ggo_visium "$@"
