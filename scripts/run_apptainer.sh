#!/usr/bin/env bash
# Deprecated wrapper — use run_container.sh instead.
# Forces Apptainer runtime explicitly.
# Usage: ./scripts/run_apptainer.sh <project_path> [command...]
SC_TOOLS_RUNTIME=apptainer exec "$(dirname "$0")/run_container.sh" "$@"
