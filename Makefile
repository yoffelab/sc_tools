# ============================================================================
# sc_tools repo-root Makefile — Developer tooling (lint, test, docs)
# ============================================================================
# Pipeline runs use Snakemake per-project (projects/<platform>/<project>/Snakefile).
# Containerized runs use Apptainer/Docker via scripts/run_container.sh.
#
# Usage:
#   make lint        - Run ruff check + format on sc_tools (required before commit)
#   make format      - Apply ruff format to sc_tools
#   make test        - Run pytest for sc_tools (sc_tools/tests/)
#   make docs        - Build Sphinx HTML docs
# ============================================================================

.PHONY: lint format test docs docs-clean docs-open wiki-sync wiki-lint

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

# Docs: build Sphinx HTML documentation (treat warnings as errors)
docs:
	sphinx-build -b html docs docs/_build/html -W --keep-going

# Docs clean: remove build artifacts and autosummary stubs
docs-clean:
	rm -rf docs/_build docs/api/generated

# Docs open: build then open in browser (macOS)
docs-open: docs
	open docs/_build/html/index.html

# Wiki: sync .gen.md files and sentinel tables from pipeline.py + registry
wiki-sync:
	python scripts/sync_wiki.py sync

# Wiki: lint .gen.md and .suggest.md files
wiki-lint:
	python scripts/sync_wiki.py lint
