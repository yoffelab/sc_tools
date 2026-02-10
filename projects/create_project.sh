#!/usr/bin/env bash
# =============================================================================
# create_project.sh - Create a new project directory under projects/<data_type>/
# =============================================================================
# Usage (from repo root):
#   ./projects/create_project.sh <project_name> <data_type>
#
# Example:
#   ./projects/create_project.sh ggo_visium visium
#   -> creates ./projects/visium/ggo_visium/ with data/, figures/, metadata/,
#      scripts/, results/, outputs/, Mission.md, Journal.md
#
# Valid data_type: visium | visium_hd | xenium | imc | cosmx
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# If script lives in projects/, repo root is parent; else script is at repo root
if [[ "$(basename "$SCRIPT_DIR")" == "projects" ]]; then
  REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
else
  REPO_ROOT="$SCRIPT_DIR"
fi
PROJECT_NAME="${1:?Usage: $0 <project_name> <data_type>}"
DATA_TYPE="${2:?Usage: $0 <project_name> <data_type>}"

VALID_TYPES=(visium visium_hd xenium imc cosmx)
if [[ ! " ${VALID_TYPES[*]} " =~ " ${DATA_TYPE} " ]]; then
  echo "Error: data_type must be one of: ${VALID_TYPES[*]}" >&2
  exit 1
fi

PROJECT_ROOT="${REPO_ROOT}/projects/${DATA_TYPE}/${PROJECT_NAME}"
SUBDIRS=(data figures metadata scripts results outputs)

mkdir -p "$PROJECT_ROOT"
for d in "${SUBDIRS[@]}"; do
  mkdir -p "${PROJECT_ROOT}/${d}"
done

# Project-specific Mission and Journal (study-specific goals and decision log)
cat > "${PROJECT_ROOT}/Mission.md" << MISSION_EOF
# Mission: ${PROJECT_NAME} (${DATA_TYPE})

**Project:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`  
**Last Updated:** (fill in)

This file holds **project-specific** goals. Repository-level pipeline and toolkit goals are in the root \`Mission.md\`.

---

## 1. Objective
(Describe the scientific objective of this project.)

---

## 2. Completed Tasks
- [ ] (Add as tasks are completed.)

---

## 3. Active Tasks (Roadmap)
(Organize by phase or theme; reference scripts under this project's \`scripts/\`.)

---

## 4. Blockers and Sanity Checks
(Project-specific blockers and statistical/data checks.)

---

## 5. Technical Decision Log (Reference this project's Journal.md)
(Key parameter and method choices for this study.)
MISSION_EOF

cat > "${PROJECT_ROOT}/Journal.md" << JOURNAL_EOF
# Research Journal & Decision Log: ${PROJECT_NAME} (${DATA_TYPE})

This journal documents the technical evolution and rationale behind **project-specific** analytical decisions. Repository-level decisions (architecture, script sanity check, scalable layout) are in the root \`Journal.md\`.

**Project:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`  
**Reference:** Root \`Mission.md\` (toolkit); this project's \`Mission.md\` (study aims).

---

## Log Entries
(Add dated entries for analysis decisions, parameter choices, and fixes specific to this project.)
JOURNAL_EOF

echo "Created project: ${PROJECT_ROOT}"
echo "  with: ${SUBDIRS[*]}, Mission.md, Journal.md"
