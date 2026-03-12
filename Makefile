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

# ============================================================================
# Container targets (Docker → Apptainer SIF → HPC deploy)
# ============================================================================
# Usage:
#   make container-build                  # Build Docker image locally
#   make container-sif                    # Convert to Apptainer SIF (Linux only)
#   make container-deploy SERVERS="brb cayuga"  # rsync SIF to HPC scratch
#   make container-all                    # Build + SIF + deploy in one step
#   make container-test                   # Quick local smoke test
# ============================================================================

DOCKER_IMAGE ?= sc_tools:latest
SIF          ?= containers/sc_tools.sif
HPC_SIF_PATH ?= /athena/elementolab/scratch/juk4007/containers/sc_tools.sif
SERVERS      ?= brb cayuga

.PHONY: container-build container-sif container-deploy container-all container-test

# Build Docker image from Dockerfile (linux/amd64 for HPC compatibility)
container-build:
	docker buildx build --platform linux/amd64 --output type=docker -t $(DOCKER_IMAGE) .
	@echo "Built $(DOCKER_IMAGE)"

# Convert Docker image to Apptainer SIF (requires apptainer; run on Linux or brb/cayuga login)
container-sif: container-build
	mkdir -p containers
	apptainer build $(SIF) docker-daemon://$(DOCKER_IMAGE)
	@echo "SIF written to $(SIF)"
	@ls -lh $(SIF)

# rsync SIF to each HPC server
container-deploy:
	@if [ ! -f $(SIF) ]; then echo "ERROR: $(SIF) not found. Run 'make container-sif' first."; exit 1; fi
	@for server in $(SERVERS); do \
	    echo "--- Deploying to $$server ---"; \
	    ssh $$server "mkdir -p $$(dirname $(HPC_SIF_PATH))"; \
	    rsync -avP --progress $(SIF) $$server:$(HPC_SIF_PATH); \
	    echo "Verifying on $$server ..."; \
	    ssh $$server "apptainer exec $(HPC_SIF_PATH) python -c 'import sc_tools, scvi, lifelines; print(\"OK on $$server\")'"; \
	done

# Full build + SIF + deploy
container-all: container-build container-sif container-deploy

# Quick local smoke test (Docker only; no HPC)
container-test:
	docker run --rm $(DOCKER_IMAGE) python -c \
	    "import sc_tools, scanpy, squidpy, scvi, lifelines, harmonypy, leidenalg; print('All imports OK')"

# Build SIF via brb (macOS workaround: apptainer is Linux-only)
# Exports as OCI tar (platform-agnostic), rsyncs to brb, builds SIF with apptainer there,
# rsyncs SIF back. OCI format avoids docker save manifest issues on macOS buildx.
OCI_TMP ?= /tmp/sc_tools_oci.tar
container-sif-via-brb: container-build
	mkdir -p containers
	@echo "--- Exporting OCI tar (linux/amd64) ---"
	docker buildx build --platform linux/amd64 --output type=oci,dest=$(OCI_TMP) .
	@echo "OCI tar size: $$(du -sh $(OCI_TMP) | cut -f1)"
	@echo "--- Uploading to brb scratch ---"
	rsync -avP --progress $(OCI_TMP) \
	    brb:/athena/elementolab/scratch/juk4007/containers/sc_tools_oci.tar
	@echo "--- Building SIF on brb ---"
	ssh brb "mkdir -p /athena/elementolab/scratch/juk4007/containers && \
	    apptainer build --force \
	        /athena/elementolab/scratch/juk4007/containers/sc_tools.sif \
	        oci-archive:///athena/elementolab/scratch/juk4007/containers/sc_tools_oci.tar && \
	    echo 'SIF built:' && ls -lh /athena/elementolab/scratch/juk4007/containers/sc_tools.sif"
	@echo "--- Rsyncing SIF back to local containers/ ---"
	rsync -avP --progress \
	    brb:/athena/elementolab/scratch/juk4007/containers/sc_tools.sif $(SIF)
	@echo "--- Cleaning up OCI tar ---"
	rm -f $(OCI_TMP)
	ssh brb "rm -f /athena/elementolab/scratch/juk4007/containers/sc_tools_oci.tar"
	@echo "Done. SIF at $(SIF)"
	@ls -lh $(SIF)

# Full macOS workflow: build Docker image, convert SIF via brb, deploy to both servers
container-all-macos: container-build container-sif-via-brb container-deploy
