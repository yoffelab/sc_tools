#!/usr/bin/env bash
# Deprecated wrapper — use run_container.sh instead.
# Forces Docker runtime explicitly.
# Usage: ./scripts/run_docker.sh <project_path> [command...]
SC_TOOLS_RUNTIME=docker exec "$(dirname "$0")/run_container.sh" "$@"
