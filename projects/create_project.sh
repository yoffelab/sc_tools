#!/usr/bin/env bash
# =============================================================================
# create_project.sh - Create a new project directory under projects/<data_type>/
# =============================================================================
# Usage (from repo root):
#   ./projects/create_project.sh <project_name> <data_type>
#
# Example:
#   ./projects/create_project.sh my_project visium
#   -> creates ./projects/visium/my_project/ with:
#      data/, figures/, metadata/, scripts/, results/, outputs/, tests/
#      Mission.md, Journal.md, journal_summary.md
#      CLAUDE.md, Snakefile, config.yaml, pyproject.toml
#
# Valid data_type: visium | visium_hd | visium_hd_cell | xenium | imc | cosmx_1k | cosmx_6k | cosmx_full_library
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(basename "$SCRIPT_DIR")" == "projects" ]]; then
  REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
else
  REPO_ROOT="$SCRIPT_DIR"
fi
PROJECT_NAME="${1:?Usage: $0 <project_name> <data_type>}"
DATA_TYPE="${2:?Usage: $0 <project_name> <data_type>}"

VALID_TYPES=(visium visium_hd visium_hd_cell xenium imc cosmx_1k cosmx_6k cosmx_full_library)
if [[ ! " ${VALID_TYPES[*]} " =~ " ${DATA_TYPE} " ]]; then
  echo "Error: data_type must be one of: ${VALID_TYPES[*]}" >&2
  exit 1
fi

PROJECT_ROOT="${REPO_ROOT}/projects/${DATA_TYPE}/${PROJECT_NAME}"
SUBDIRS=(data figures metadata scripts results outputs tests)

mkdir -p "$PROJECT_ROOT"
for d in "${SUBDIRS[@]}"; do
  mkdir -p "${PROJECT_ROOT}/${d}"
done

# ---- Mission.md ----
cat > "${PROJECT_ROOT}/Mission.md" << MISSION_EOF
# Mission: ${PROJECT_NAME} (${DATA_TYPE})

**Project:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`
**Current Status:** Setup
**Last Updated:** $(date +%Y-%m-%d)

Project-specific goals. Repository-level pipeline and toolkit goals are in the root \`Mission.md\`.

---

## 1. Objective
(Describe the scientific objective of this project.)

---

## 2. Phase Alignment

| Phase | Status | Key Tasks |
|-------|--------|-----------|
| **1** | Pending | Data ingestion, QC |
| **2** | Pending | Metadata attachment |
| **3** | Pending | Preprocessing, clustering |
| **3.5** | Pending | Demographics |
| **3.5b** | Pending | Gene scoring, cell typing, deconvolution |
| **4** | Pending | Manual cell typing |
| **5** | Pending | Downstream biology |
| **6-7** | Pending | Meta analysis |

---

## 3. Active Tasks

- [ ] Phase 1: data ingestion

---

## 4. Blockers and Sanity Checks
(Project-specific blockers and checks.)

---

## 5. Technical Decisions
(Key parameter and method choices — reference Journal.md for full log.)
MISSION_EOF

# ---- Journal.md ----
cat > "${PROJECT_ROOT}/Journal.md" << JOURNAL_EOF
# Research Journal: ${PROJECT_NAME} (${DATA_TYPE})

**Project:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`

---

## Log Entries

### $(date +%Y-%m-%d) — Project created
- Scaffolded with create_project.sh.
JOURNAL_EOF

# ---- journal_summary.md ----
cat > "${PROJECT_ROOT}/journal_summary.md" << SUMMARY_EOF
# Journal Summary: ${PROJECT_NAME} (${DATA_TYPE})

Condensed summary of \`Journal.md\`. Full entries in Journal.md.

## Project scope
(One short paragraph on scientific goal and current phase.)

## Recent phase
(Bullets on latest decisions and outcomes.)

## Key conventions
(One line each: panel names, checkpoint names, etc.)
SUMMARY_EOF

# ---- CLAUDE.md ----
cat > "${PROJECT_ROOT}/CLAUDE.md" << CLAUDE_EOF
# ${PROJECT_NAME} — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @journal_summary.md — recent decisions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Platform:** ${DATA_TYPE}
**Path:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`

(Describe scientific objective here.)

---

## Running This Project

\`\`\`bash
# From repo root:
./scripts/run_container.sh projects/${DATA_TYPE}/${PROJECT_NAME} python scripts/<script>.py
snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile <target>

# Local conda (no container):
conda activate ${PROJECT_NAME}
SC_TOOLS_RUNTIME=none snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile <target>

# Tests:
pytest projects/${DATA_TYPE}/${PROJECT_NAME}/tests/ -v

# Create conda env (one-time setup, from repo root):
#   conda create -n ${PROJECT_NAME} python=3.10 -y && conda activate ${PROJECT_NAME}
#   uv pip install -e ".[deconvolution]"
\`\`\`

---

## Key Files

| Path | Description |
|------|-------------|
| \`results/adata.raw.p1.h5ad\` | Phase 1 output |
| \`results/adata.normalized.scored.p35.h5ad\` | Phase 3.5b output (primary analysis input) |
| \`metadata/gene_signatures.json\` | Gene signatures |
| \`figures/manuscript/\` | Publication figures |
CLAUDE_EOF

# ---- config.yaml ----
cat > "${PROJECT_ROOT}/config.yaml" << CONFIG_EOF
# ${PROJECT_NAME} Snakemake config
repo_root: "../../.."
project_rel: "projects/${DATA_TYPE}/${PROJECT_NAME}"
container_sif: "containers/sc_tools.sif"
CONFIG_EOF

# ---- Snakefile ----
cat > "${PROJECT_ROOT}/Snakefile" << SNAKEFILE_EOF
# ${PROJECT_NAME} (${DATA_TYPE}) Snakemake pipeline
# Run from repo root: snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile [target]
# Or from project dir: snakemake -d . -s Snakefile [target]
# Runtime: auto-detected via run_container.sh (Docker on macOS, Apptainer on Linux).
# Local conda: conda activate ${PROJECT_NAME} && SC_TOOLS_RUNTIME=none snakemake -d . -s Snakefile <target>

configfile: "config.yaml"

import os
import pandas as pd

ROOT = config["repo_root"]
PROJECT = config["project_rel"]

def run_container(script, args=""):
    cmd = f'cd {ROOT} && ./scripts/run_container.sh {PROJECT} python {script}'
    return cmd + (f' {args}' if args else '')

# ---- Phase 0: Upstream Raw Data Processing ----
PHASE0_DIR = "metadata/phase0"

def load_phase0_manifest():
    all_tsv = f"{PHASE0_DIR}/all_samples.tsv"
    if os.path.exists(all_tsv):
        return pd.read_csv(all_tsv, sep="\t")
    return pd.DataFrame()

PHASE0_MANIFEST = load_phase0_manifest()
PHASE0_SAMPLES = list(PHASE0_MANIFEST["sample_id"]) if len(PHASE0_MANIFEST) > 0 else []

rule phase0:
    input: expand("data/{sample_id}/outs", sample_id=PHASE0_SAMPLES)

# ---- Phase 1: Data Ingestion & QC ----
rule adata_raw_p1:
    output: "results/adata.raw.p1.h5ad"
    input: "scripts/ingest.py"
    shell: run_container("scripts/ingest.py")

rule validate_p1:
    input: "results/adata.raw.p1.h5ad"
    output: touch("results/.adata.raw.p1.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.raw.p1.h5ad --phase p1 --fix"

rule phase1:
    input: "results/adata.raw.p1.h5ad", "results/.adata.raw.p1.validated"

# ---- Phase 2: Metadata Attachment ----
rule adata_annotated_p2:
    output: "results/adata.annotated.p2.h5ad"
    input: "results/adata.raw.p1.h5ad", "results/.adata.raw.p1.validated", "metadata/sample_metadata.csv"
    shell: run_container("scripts/attach_metadata.py")

rule validate_p2:
    input: "results/adata.annotated.p2.h5ad"
    output: touch("results/.adata.annotated.p2.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.annotated.p2.h5ad --phase p2 --fix"

rule phase2:
    input: "results/adata.annotated.p2.h5ad", "results/.adata.annotated.p2.validated"

# ---- Phase 3: Preprocessing ----
rule adata_normalized_p3:
    output: "results/adata.normalized.p3.h5ad"
    input: "results/adata.annotated.p2.h5ad", "results/.adata.annotated.p2.validated", "scripts/preprocess.py"
    shell: run_container("scripts/preprocess.py")

rule validate_p3:
    input: "results/adata.normalized.p3.h5ad"
    output: touch("results/.adata.normalized.p3.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.normalized.p3.h5ad --phase p3 --fix"

rule phase3:
    input: "results/adata.normalized.p3.h5ad", "results/.adata.normalized.p3.validated"

# ---- Phase 3.5b: Gene Scoring ----
rule adata_p35:
    output: "results/adata.normalized.scored.p35.h5ad"
    input: "results/adata.normalized.p3.h5ad", "results/.adata.normalized.p3.validated", "metadata/gene_signatures.json", "scripts/score_signatures.py"
    shell: run_container("scripts/score_signatures.py")

rule validate_p35:
    input: "results/adata.normalized.scored.p35.h5ad"
    output: touch("results/.adata.normalized.scored.p35.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.normalized.scored.p35.h5ad --phase p35 --fix"

rule phase35b:
    input: "results/adata.normalized.scored.p35.h5ad", "results/.adata.normalized.scored.p35.validated"

# ---- Phase 4: Cell Typing ----
rule adata_celltyped_p4:
    output: "results/adata.celltyped.p4.h5ad"
    input: "results/adata.normalized.scored.p35.h5ad", "results/.adata.normalized.scored.p35.validated", "metadata/celltype_map.json", "scripts/apply_celltype.py"
    shell: run_container("scripts/apply_celltype.py")

rule validate_p4:
    input: "results/adata.celltyped.p4.h5ad"
    output: touch("results/.adata.celltyped.p4.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.celltyped.p4.h5ad --phase p4 --fix"

rule phase4:
    input: "results/adata.celltyped.p4.h5ad", "results/.adata.celltyped.p4.validated"

# ---- QC reports (date-versioned HTML) ----
rule qc_pre_filter:
    input: "results/adata.raw.p1.h5ad"
    output: touch("figures/QC/pre_filter_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report pre_filter --adata results/adata.raw.p1.h5ad --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_post_filter:
    input: "results/adata.raw.p1.h5ad", "results/adata.annotated.p2.h5ad"
    output: touch("figures/QC/post_filter_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report post_filter --adata results/adata.raw.p1.h5ad"
        + " --adata-post results/adata.annotated.p2.h5ad --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_post_integration:
    input: "results/adata.normalized.p3.h5ad"
    output: touch("figures/QC/post_integration_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report post_integration --adata-integrated results/adata.normalized.p3.h5ad"
        + " --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_report:
    input: rules.qc_pre_filter.output, rules.qc_post_filter.output, rules.qc_post_integration.output
    output: touch("figures/QC/qc_report.done")

# ---- Default ----
rule all:
    input: "results/adata.normalized.scored.p35.h5ad"
SNAKEFILE_EOF

# ---- pyproject.toml (per-project package descriptor) ----
cat > "${PROJECT_ROOT}/pyproject.toml" << PYPROJECT_EOF
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"

[project]
name = "${PROJECT_NAME}"
version = "0.1.0"
description = "${PROJECT_NAME} (${DATA_TYPE}) analysis (scripts only)"
requires-python = ">=3.10"
# Add project-specific deps here (beyond sc-tools base + deconvolution extra).
# Most deps come from root: uv pip install -e ".[deconvolution]"
dependencies = []

[tool.setuptools]
packages = []  # scripts-only project; no Python packages to install
PYPROJECT_EOF

echo "Created project: ${PROJECT_ROOT}"
echo "  dirs:  ${SUBDIRS[*]}"
echo "  files: Mission.md, Journal.md, journal_summary.md, CLAUDE.md"
echo "         Snakefile, config.yaml, pyproject.toml"
echo ""
echo "Next steps:"
echo "  1. Edit Mission.md with your scientific objective"
echo "  2. Add ingestion script: ${PROJECT_ROOT}/scripts/ingest.py"
echo "  3. Create conda env:"
echo "       conda create -n ${PROJECT_NAME} python=3.10 -y"
echo "       conda activate ${PROJECT_NAME}"
echo "       uv pip install -e '.[deconvolution]'   # from repo root"
echo "  4. Run (container): snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile phase1"
echo "  5. Run (local):     conda activate ${PROJECT_NAME} && SC_TOOLS_RUNTIME=none snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile phase1"
