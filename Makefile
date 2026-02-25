# ============================================================================
# sc_tools repo-root Makefile — Lint, test, and other repo-wide targets
# ============================================================================
# DEPRECATED: Prefer Snakemake (Snakefile at repo root: snakemake lint, test, all).
# Containerized runs use Apptainer; see project_setup.md.
#
# Pipeline Makefile lives in projects/visium/ggo_visium/Makefile (run with
# make -f projects/visium/ggo_visium/Makefile or PROJECT=...).
#
# Usage:
#   make lint        - Run ruff check + format on sc_tools (required before commit)
#   make format      - Apply ruff format to sc_tools
#   make test        - Run pytest for sc_tools (sc_tools/tests/)
# ============================================================================

.PHONY: lint format test

# Lint: check sc_tools (package). Scripts: add once fixed.
lint:
	ruff check sc_tools
	ruff format --check sc_tools

# Format: apply ruff format to sc_tools
format:
	ruff format sc_tools

# Test: run pytest for sc_tools package (from repo root so project paths resolve)
# Uses python -m pytest so the active Python env (conda/venv) is used
test:
	python -m pytest sc_tools/tests -v
