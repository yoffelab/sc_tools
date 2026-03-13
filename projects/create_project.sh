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
#      figures/, metadata/, scripts/, results/, outputs/, tests/
#      (data/ is NOT scaffolded — data paths live in registry.db projects table)
#      Plan.md, Journal.md
#      Snakefile, config.yaml, pyproject.toml
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
SUBDIRS=(metadata scripts results outputs tests)

mkdir -p "$PROJECT_ROOT"
for d in "${SUBDIRS[@]}"; do
  mkdir -p "${PROJECT_ROOT}/${d}"
done

# figures/ subdirectory layout
# scratch/       — throwaway exploratory plots; never evaluated by the hook
# QC/            — auto-generated pipeline QC reports; technical checks only
# exploratory/   — analysis in progress; full eval, standard bar
# insightful/    — clear biological findings, pre-manuscript; full eval, high bar
# manuscript/    — final publication figures; full eval, strictest bar
# supplementary/ — supplementary material; same bar as manuscript
for fig_subdir in scratch QC exploratory insightful manuscript supplementary; do
  mkdir -p "${PROJECT_ROOT}/figures/${fig_subdir}"
done

# ---- Phase 0 status template ----
mkdir -p "${PROJECT_ROOT}/metadata/phase0"
printf "sample_id\tstatus\tnotes\n" > "${PROJECT_ROOT}/metadata/phase0/phase0_status.tsv"

# ---- Plan.md ----
cat > "${PROJECT_ROOT}/Plan.md" << PLAN_EOF
# Plan: ${PROJECT_NAME} (${DATA_TYPE})

## Context

(Scientific objective of this project.)

**Project:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`
**Platform:** ${DATA_TYPE}

---

## Phase Status

Query registry via \`get_phase_status\` or check registry. Last completed: (none).

---

## Tasks

### Current Priority

- [ ] (First task)

### Backlog

- [ ] (Future tasks)

---

## Blocked

(Items waiting on external input, data, or decisions.)

---

## Technical Decisions

(Key parameter choices, method selections, conventions specific to this project.
Reference Journal.md for the full dated log.)

---

## Key Files

| Path | Description |
|------|-------------|
| \`results/adata.filtered.h5ad\` | qc_filter output |
| \`results/adata.annotated.h5ad\` | metadata_attach output |
| \`results/adata.normalized.h5ad\` | preprocess output |
| \`results/adata.scored.h5ad\` | scoring output |
| \`results/adata.celltyped.h5ad\` | celltype_manual output |
| \`metadata/gene_signatures.json\` | Gene signatures |

---

## Conventions

(Project-specific conventions: tumor stage ordering, comparison mode,
normalization parameters, panel names, etc. Repo-wide conventions are
in root CLAUDE.md and docs/skills.md -- do not duplicate here.)
PLAN_EOF

# ---- Journal.md ----
cat > "${PROJECT_ROOT}/Journal.md" << JOURNAL_EOF
# Research Journal: ${PROJECT_NAME} (${DATA_TYPE})

**Project:** \`projects/${DATA_TYPE}/${PROJECT_NAME}\`

---

## Log Entries

### $(date +%Y-%m-%d) — Project created
- Scaffolded with create_project.sh.
JOURNAL_EOF

# ---- config.yaml ----
cat > "${PROJECT_ROOT}/config.yaml" << CONFIG_EOF
# ${PROJECT_NAME} Snakemake config
repo_root: "../../.."
project_rel: "projects/${DATA_TYPE}/${PROJECT_NAME}"
container_sif: "containers/sc_tools.sif"

# Phase 0: upstream tool paths (adjust for your HPC)
# spaceranger_path: "/path/to/spaceranger"
# transcriptome: "/path/to/refdata-gex-GRCh38-2024-A"
# probe_set: "/path/to/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv"
# output_dir: "spaceranger_outputs"
# create_bam: false
# localcores: 32
# localmem: 220

# SLURM resource defaults (used by generate_phase0_sbatch.py)
# slurm:
#   partition: "scu-cpu"
#   cpus_per_task: 32
#   mem: "240G"
#   time: "2-00:00:00"
CONFIG_EOF

# ---- Snakefile ----
cat > "${PROJECT_ROOT}/Snakefile" << SNAKEFILE_EOF
# ${PROJECT_NAME} (${DATA_TYPE}) Snakemake pipeline
# Run from repo root: snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile [target]
# Or from project dir: snakemake -d . -s Snakefile [target]
# Runtime: auto-detected via run_container.sh (Docker on macOS, Apptainer on Linux).
# Local conda: conda activate ${PROJECT_NAME} && SC_TOOLS_RUNTIME=none snakemake -d . -s Snakefile <target>
#
# ── Phase naming convention ───────────────────────────────────────────────────
# This Snakefile uses SEMANTIC SLUG names (new nomenclature) defined in
# sc_tools.pipeline. Rule names and checkpoint filenames follow the slugs:
#
#   Slug              Old code (DEPRECATED)   Checkpoint file
#   qc_filter         p1                      results/adata.filtered.h5ad
#   metadata_attach   p2                      results/adata.annotated.h5ad
#   preprocess        p3                      results/adata.normalized.h5ad
#   scoring           p35                     results/adata.scored.h5ad
#   celltype_manual   p4                      results/adata.celltyped.h5ad
#
# The old p1/p2/p3/p35/p4 names and adata.*.p1.h5ad filenames are OLD
# NOMENCLATURE. They are still accepted by sc_tools.validate for migration
# purposes but will emit DeprecationWarning. New projects should use slugs.
# ─────────────────────────────────────────────────────────────────────────────

configfile: "config.yaml"

import os
import pandas as pd

ROOT = config["repo_root"]
PROJECT = config["project_rel"]

def run_container(script, args=""):
    cmd = f'cd {ROOT} && ./scripts/run_container.sh {PROJECT} python {script}'
    return cmd + (f' {args}' if args else '')

# ---- ingest_raw / ingest_load: Upstream Raw Data Processing (Phase 0a/0b) ----
PHASE0_DIR = "metadata/phase0"

def load_phase0_manifest():
    all_tsv = f"{PHASE0_DIR}/all_samples.tsv"
    if os.path.exists(all_tsv):
        return pd.read_csv(all_tsv, sep="\t")
    return pd.DataFrame()

PHASE0_MANIFEST = load_phase0_manifest()
PHASE0_SAMPLES = list(PHASE0_MANIFEST["sample_id"]) if len(PHASE0_MANIFEST) > 0 else []

# Data paths are resolved from registry.db (not a local data/ directory).
# Override DATA_DIR in config.yaml or query via: sc_tools.registry.get_project(project).data_dir
rule ingest_raw:
    input: expand(config.get("data_dir", "data") + "/{sample_id}/outs", sample_id=PHASE0_SAMPLES)

# ---- qc_filter: QC + Concatenation ----
rule adata_qc_filter:
    output: "results/adata.filtered.h5ad"
    input: "scripts/ingest.py"
    shell: run_container("scripts/ingest.py")

rule validate_qc_filter:
    input: "results/adata.filtered.h5ad"
    output: touch("results/.adata.filtered.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.filtered.h5ad --phase qc_filter --fix"

rule qc_filter:
    input: "results/adata.filtered.h5ad", "results/.adata.filtered.validated"

# ---- metadata_attach: Metadata Attachment ----
rule adata_metadata_attach:
    output: "results/adata.annotated.h5ad"
    input: "results/adata.filtered.h5ad", "results/.adata.filtered.validated", "metadata/sample_metadata.csv"
    shell: run_container("scripts/attach_metadata.py")

rule validate_metadata_attach:
    input: "results/adata.annotated.h5ad"
    output: touch("results/.adata.annotated.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.annotated.h5ad --phase metadata_attach --fix"

rule metadata_attach:
    input: "results/adata.annotated.h5ad", "results/.adata.annotated.validated"

# ---- preprocess: Normalize + Integrate + Cluster ----
rule adata_preprocess:
    output: "results/adata.normalized.h5ad"
    input: "results/adata.annotated.h5ad", "results/.adata.annotated.validated", "scripts/preprocess.py"
    shell: run_container("scripts/preprocess.py")

rule validate_preprocess:
    input: "results/adata.normalized.h5ad"
    output: touch("results/.adata.normalized.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.normalized.h5ad --phase preprocess --fix"

rule preprocess:
    input: "results/adata.normalized.h5ad", "results/.adata.normalized.validated"

# ---- scoring: Gene Scoring + Auto Cell Typing ----
rule adata_scoring:
    output: "results/adata.scored.h5ad"
    input: "results/adata.normalized.h5ad", "results/.adata.normalized.validated", "metadata/gene_signatures.json", "scripts/score_signatures.py"
    shell: run_container("scripts/score_signatures.py")

rule validate_scoring:
    input: "results/adata.scored.h5ad"
    output: touch("results/.adata.scored.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.scored.h5ad --phase scoring --fix"

rule scoring:
    input: "results/adata.scored.h5ad", "results/.adata.scored.validated"

# ---- celltype_manual: Manual Cell Typing (iterative, optional) ----
rule adata_celltype_manual:
    output: "results/adata.celltyped.h5ad"
    input: "results/adata.scored.h5ad", "results/.adata.scored.validated", "metadata/celltype_map.json", "scripts/apply_celltype.py"
    shell: run_container("scripts/apply_celltype.py")

rule validate_celltype_manual:
    input: "results/adata.celltyped.h5ad"
    output: touch("results/.adata.celltyped.validated")
    shell: f"cd {{ROOT}} && python scripts/validate_checkpoint.py {{PROJECT}}/results/adata.celltyped.h5ad --phase celltype_manual --fix"

rule celltype_manual:
    input: "results/adata.celltyped.h5ad", "results/.adata.celltyped.validated"

# ---- QC reports (date-versioned HTML) ----
rule qc_pre_filter:
    input: "results/adata.filtered.h5ad"
    output: touch("figures/QC/pre_filter_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report pre_filter --adata results/adata.filtered.h5ad --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_post_filter:
    input: "results/adata.filtered.h5ad", "results/adata.annotated.h5ad"
    output: touch("figures/QC/post_filter_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report post_filter --adata results/adata.filtered.h5ad"
        + " --adata-post results/adata.annotated.h5ad --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_post_integration:
    input: "results/adata.normalized.h5ad"
    output: touch("figures/QC/post_integration_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report post_integration --adata-integrated results/adata.normalized.h5ad"
        + " --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_post_celltyping:
    input: "results/adata.celltyped.h5ad", "results/.adata.celltyped.validated"
    output: touch("figures/QC/post_celltyping_qc.done")
    shell: (
        run_container(ROOT + "/scripts/run_qc_report.py")
        + " --report post_celltyping --adata-integrated results/adata.celltyped.h5ad"
        + " --figures-dir figures --modality ${DATA_TYPE}"
    )

rule qc_report:
    input: rules.qc_pre_filter.output, rules.qc_post_filter.output, rules.qc_post_integration.output
    output: touch("figures/QC/qc_report.done")

# ---- Default ----
rule all:
    input: "results/adata.scored.h5ad"
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
echo "  dirs:  ${SUBDIRS[*]} figures/{scratch,QC,exploratory,insightful,manuscript,supplementary}"
echo "  files: Plan.md, Journal.md"
echo "         Snakefile, config.yaml, pyproject.toml"
echo ""
echo "Next steps:"
echo "  1. Edit Plan.md with your scientific objective"
echo "  2. Add ingestion script: ${PROJECT_ROOT}/scripts/ingest.py"
echo "  3. Create conda env:"
echo "       conda create -n ${PROJECT_NAME} python=3.11 -y"
echo "       conda activate ${PROJECT_NAME}"
echo "       uv pip install -e '.[deconvolution]'   # from repo root"
echo "  4. Run (container): snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile qc_filter"
echo "  5. Run (local):     conda activate ${PROJECT_NAME} && SC_TOOLS_RUNTIME=none snakemake -d projects/${DATA_TYPE}/${PROJECT_NAME} -s projects/${DATA_TYPE}/${PROJECT_NAME}/Snakefile qc_filter"
