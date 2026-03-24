# Roadmap: sc_tools Agent-Native CLI

## Overview

Transform sc_tools from a library-with-scripts into an agent-native CLI (`sct`) with structured JSON output, file-based provenance, and memory-safe data access. The journey starts by fixing known internals bugs (immediate value for active projects), then builds the CLI foundation and commands layer-by-layer, and finishes with scientific analysis gaps and multi-omic assembly. Phases 1 and 2 can run in parallel since they have no mutual dependencies.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 1: Benchmark Fixes** - Fix integration benchmark bugs and add missing features for immediate project value
- [x] **Phase 2: CLI Foundation** - Register `sct` entry point with Typer, CLIResult envelope, structured output, and error handling (completed 2026-03-21)
- [ ] **Phase 3: Core Commands** - Wrap existing sc_tools operations as CLI commands with shared Result type
- [ ] **Phase 4: CLI Discovery** - Self-describing command catalog and schema export for agent introspection
- [ ] **Phase 5: Provenance** - File-based provenance tracking with JSON sidecars and lineage queries
- [ ] **Phase 6: Scientific Gaps** - Pseudobulk DE, marker validation, subject-level metadata, panel-aware cell typing
- [ ] **Phase 7: Memory Safety** - Tiered IO gateway, memory estimation, and dry-run mode for large datasets
- [ ] **Phase 8: Multi-Omic Assembly** - Cross-modality metadata join, MuData assembly, and joint embedding

## Phase Details

### Phase 1: Benchmark Fixes
**Goal**: Integration benchmarking produces correct, reproducible, memory-efficient results on real project data
**Depends on**: Nothing (can start immediately)
**Requirements**: BM-01, BM-02, BM-03, BM-04, BM-05, BM-06, BM-07, TST-01, TST-02, TST-03
**Success Criteria** (what must be TRUE):
  1. `compare_integrations()` loads pre-computed h5ad embeddings via h5py without loading full AnnData (peak memory <2GB for 2.5M-cell dataset)
  2. Running benchmarks on embeddings with NaN rows (e.g., resolVI) produces valid metrics without crashing or silent NaN propagation
  3. `_stratified_subsample` preserves group proportions (no truncation bias toward low-index cells)
  4. Benchmark output includes `runtime_s` column and all parameter provenance (batch_weight, bio_weight, seed, resolution)
  5. `_recipe_targeted_panel` skips normalization when scVI integration is selected (raw counts preserved for scVI)
**Plans**: 2 plans

Plans:
- [x] 01-01-PLAN.md — Fix _stratified_subsample, NaN masking, subsample_n, runtime_s, param provenance (BM-02..06, TST-01..03)
- [x] 01-02-PLAN.md — h5py embedding loading, bio_key parameter, targeted panel scVI fix (BM-01, BM-07)

### Phase 2: CLI Foundation
**Goal**: `sct` is an installable, fast-starting CLI that produces structured JSON output and handles errors with actionable taxonomy
**Depends on**: Nothing (can run in parallel with Phase 1)
**Requirements**: CLI-01, CLI-02, CLI-03, CLI-04, CLI-05, CLI-06, CLI-07, CLI-08, TST-04
**Success Criteria** (what must be TRUE):
  1. `pip install -e .` registers `sct` and `which sct` returns a path
  2. `sct help` returns in under 500ms (no heavy imports at startup)
  3. Every command outputs a valid JSON CLIResult envelope to stdout with status, command, data, artifacts, provenance, and message fields
  4. `sct --human <cmd>` renders Rich-formatted output to stderr while keeping JSON on stdout
  5. Commands exit with semantic codes (0/1/2/3) and errors include retryable/fixable/fatal classification with actionable suggestions
**Plans**: 2 plans

Plans:
- [x] 02-01-PLAN.md — CLIResult Pydantic model, exception hierarchy, pyproject.toml cli deps (CLI-03, CLI-05, CLI-06)
- [x] 02-02-PLAN.md — Typer app, stub command groups, error handler, --human, tests, MCP proof-of-concept (CLI-01, CLI-02, CLI-04, CLI-07, CLI-08, TST-04)

### Phase 3: Core Commands
**Goal**: Agents can run QC, preprocessing, validation, benchmarking, status checks, and report generation through `sct` commands instead of ad-hoc scripts
**Depends on**: Phase 2 (needs CLIResult envelope and output machinery)
**Requirements**: CMD-01, CMD-02, CMD-03, CMD-04, CMD-05, CMD-06, CMD-07, CMD-08, TST-05, TST-06
**Success Criteria** (what must be TRUE):
  1. `sct qc run <file>` produces JSON metrics summary; `sct preprocess run <file>` dispatches modality-aware recipe and writes output h5ad
  2. `sct benchmark integration --from-dir <dir>` loads pre-computed embeddings, computes comparison metrics, and outputs ranked JSON
  3. `sct validate <phase> <file>` checks checkpoint against PhaseSpec and returns pass/fail with specific failures listed
  4. `sct status` shows pipeline DAG state (completed phases, available next, checkpoint paths) from registry
  5. CLI commands and MCP tools share the same Result type (single implementation, dual serialization verified by test)
**Plans**: 4 plans

Plans:
- [x] 03-01-PLAN.md — Migrate cli.py to cli/ package, add _check_deps utility, create test fixtures (CMD-08, TST-05)
- [x] 03-02-PLAN.md — Implement validate, status, qc run, and report commands (CMD-01, CMD-03, CMD-05, CMD-06)
- [x] 03-03-PLAN.md — Implement preprocess run and benchmark integration commands (CMD-02, CMD-04)
- [x] 03-04-PLAN.md — CLI integration tests, E2E test scaffold, shared Result verification (CMD-07, TST-05, TST-06)

### Phase 4: CLI Discovery
**Goal**: Agents can programmatically discover all available commands, their parameters, and output schemas without parsing help text
**Depends on**: Phase 3 (needs commands to introspect)
**Requirements**: DSC-01, DSC-02, DSC-03
**Success Criteria** (what must be TRUE):
  1. `sct list-commands --json` returns a machine-readable catalog with every command's name, params, types, and defaults
  2. `sct describe <cmd>` returns JSON schema for a specific command's input params and output format
  3. `sct schema` exports the full CLI contract as a single JSON document (command tree + Pydantic output schemas)
**Plans**: 1 plan

Plans:
- [ ] 04-01-PLAN.md — Discovery commands: list-commands, describe, schema with Typer/Click introspection (DSC-01, DSC-02, DSC-03)

### Phase 5: Provenance
**Goal**: Every CLI output has a traceable lineage back to its inputs, parameters, and environment
**Depends on**: Phase 2 (needs CLI output machinery for sidecar writing)
**Requirements**: PRV-01, PRV-02, PRV-03, PRV-04, PRV-05
**Success Criteria** (what must be TRUE):
  1. Every data-producing CLI command writes a `.provenance.json` sidecar alongside output files
  2. Sidecar contains command, params, input files with SHA256 checksums, sc_tools version, timestamp, runtime_s, and peak_memory_mb
  3. `sct provenance show <file>` displays the provenance record; `sct provenance trace <file>` walks input references to build a full lineage DAG
  4. Leiden clustering results include resolution and random_state in provenance, and identical parameters produce identical results
**Plans**: 2 plans

Plans:
- [x] 05-01-PLAN.md — Provenance model, utilities, sidecar hook in cli_handler, Leiden random_state (PRV-01, PRV-02, PRV-05)
- [ ] 05-02-PLAN.md — Provenance show and trace CLI commands, lineage trace engine (PRV-03, PRV-04)

### Phase 6: Scientific Gaps
**Goal**: sc_tools supports pseudobulk DE, marker validation, subject-level metadata, and panel-aware cell typing for rigorous multi-sample analysis
**Depends on**: Phase 1 (needs fixed internals for stable wrapping)
**Requirements**: SCI-01, SCI-02, SCI-03, SCI-04
**Success Criteria** (what must be TRUE):
  1. `sc_tools.tl.de` aggregates counts by subject_id + celltype and runs PyDESeq2 with batch covariates, producing a per-gene results table
  2. After cell typing, a marker validation report (dotplot/heatmap) shows top markers per assigned type and flags types with low canonical marker expression
  3. `subject_id` is enforced as distinct from `library_id` at ingestion for multi-sample projects, and batch-condition confounding is validated at registration
  4. When `n_vars < 1000` (targeted panels), cell typing dispatch restricts to panel-validated methods and warns if a whole-transcriptome model is applied
**Plans**: 3 plans

Plans:
- [x] 06-01-PLAN.md — Subject metadata validation and panel-aware cell typing dispatch (SCI-03, SCI-04)
- [x] 06-02-PLAN.md — Pseudobulk DE module with PyDESeq2 and sct de run CLI command (SCI-01)
- [ ] 06-03-PLAN.md — Marker validation report integrated into post-celltyping QC (SCI-02)

### Phase 7: Memory Safety
**Goal**: Large datasets (2.5M cells, 25G h5ad) can be processed without OOM through tiered loading, pre-execution estimation, and dry-run validation
**Depends on**: Phase 2 (needs CLI infrastructure for --dry-run flag and estimate command)
**Requirements**: MEM-01, MEM-02, MEM-03
**Success Criteria** (what must be TRUE):
  1. IO Gateway loads metadata via h5py, summaries via backed mode, and full data only for compute -- peak memory for metadata queries on 25G files stays under 1GB
  2. `sct estimate <command> <args>` returns projected peak memory and runtime based on cell/gene count and method before execution starts
  3. `--dry-run` on any data-touching command validates inputs, reports planned operations, and exits without modifying data
**Plans**: TBD

Plans:
- [ ] 07-01: TBD
- [ ] 07-02: TBD

### Phase 8: Multi-Omic Assembly
**Goal**: Independently processed modalities (scRNA, IMC, Visium, Xenium) can be assembled into a unified MuData object for cross-modal patient-level analysis
**Depends on**: Phase 6 (needs subject_id metadata model)
**Requirements**: MOM-01, MOM-02, MOM-03, MOM-04
**Success Criteria** (what must be TRUE):
  1. Patient/subject metadata joins across modalities using validated subject_id, producing a unified metadata table
  2. MuData assembly combines independently processed per-modality AnnData objects into a single MuData with shared obs
  3. Joint embedding via MOFA+, multiVI, or totalVI produces a shared latent space across modalities
  4. Cross-modal queries work at patient/project level (e.g., cell type proportions for a patient across all modalities)
**Plans**: TBD

Plans:
- [ ] 08-01: TBD
- [ ] 08-02: TBD

## Progress

**Execution Order:**
Phases execute in numeric order. Phases 1 and 2 have no mutual dependencies and can be planned/executed in parallel.

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Benchmark Fixes | 0/2 | Planning complete | - |
| 2. CLI Foundation | 2/2 | Complete   | 2026-03-21 |
| 3. Core Commands | 0/4 | Planning complete | - |
| 4. CLI Discovery | 0/1 | Planning complete | - |
| 5. Provenance | 1/2 | In Progress|  |
| 6. Scientific Gaps | 2/3 | In Progress|  |
| 7. Memory Safety | 0/2 | Not started | - |
| 8. Multi-Omic Assembly | 0/2 | Not started | - |
