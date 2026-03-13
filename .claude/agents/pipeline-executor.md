---
name: pipeline-executor
description: Run a pipeline phase end-to-end — load data, execute, validate checkpoint, report.
skills: [sc-tools-skills, verification-before-completion]
tools_expected: [Read, Write, Edit, Bash, Glob, Grep]
---

# Pipeline Executor Agent

Runs a single pipeline phase from input checkpoint to output checkpoint.

## Required context in brief
- Project path: `projects/<platform>/<project>/`
- Phase to execute (e.g., `qc_filter`, `preprocess`, `scoring`)
- Input checkpoint path and confirmation it exists
- Expected output checkpoint path
- Validation command: `python scripts/validate_checkpoint.py <output> --phase <phase>`

## Standards to inline
From `docs/Architecture.md` §2.2 — the checkpoint metadata requirements for the target phase:
- What obs/obsm/uns keys must exist
- What X should contain (raw counts vs normalized)
- Whether samples should be concatenated

From `docs/skills.md` (via sc-tools-skills):
- FDR correction: Benjamini-Hochberg always
- Signature scores in obsm, not obs
- Colors in uns; never overwrite if length matches

## How to run scripts
- Local: `./scripts/run_container.sh projects/<platform>/<project> python scripts/foo.py`
- Snakemake: `snakemake -d projects/<platform>/<project> -s projects/<platform>/<project>/Snakefile <target>`
- Full container/runtime details: `docs/project_setup.md`

## Verification
Agent must run `validate_checkpoint.py` on output before reporting success.
Report: n_obs, n_vars, n_samples, any warnings.

Repo root: see docs/Architecture.md section 1 for directory layout.
