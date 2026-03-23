# Requirements: sc_tools Agent-Native CLI

**Defined:** 2026-03-20
**Core Value:** Agents never write throwaway scripts — every comp bio operation is callable via a stable CLI with structured I/O

## v1 Requirements

Requirements for initial release. Each maps to roadmap phases.

### Benchmark Fixes

- [x] **BM-01**: Load pre-computed integration embeddings from separate per-method h5ad files using h5py (no full AnnData load; memory <2GB for 2.5M cells)
- [x] **BM-02**: Filter NaN rows per-embedding before metric computation in `compare_integrations()` (resolVI produces NaN for cells with <5 HVG counts)
- [x] **BM-03**: Configurable subsampling via `subsample_n` parameter in `compare_integrations()` (default 50K, using `_stratified_subsample`)
- [x] **BM-04**: Fix `_stratified_subsample` truncation bias — replace `sorted(indices)[:n]` with proportional group downsampling
- [x] **BM-05**: Add `runtime_s` column to `run_integration_benchmark` output for method comparison
- [x] **BM-06**: Store benchmark parameters (batch_weight, bio_weight, seed, resolution, scib backend) alongside results
- [x] **BM-07**: Fix `_recipe_targeted_panel` — skip normalization when scVI integration is selected (currently normalizes then warns about raw counts)

### CLI Foundation

- [x] **CLI-01**: Register `sct` entry point in `pyproject.toml` `[project.scripts]` section
- [x] **CLI-02**: Typer-based CLI app with command groups mapping to pipeline phases (qc, preprocess, integrate, benchmark, celltype)
- [x] **CLI-03**: Pydantic CLIResult envelope for all commands: `{status, command, data, artifacts, provenance, message}`
- [x] **CLI-04**: JSON output to stdout by default; `--human` flag renders Rich-formatted tables/text to stderr
- [x] **CLI-05**: Semantic exit codes: 0=success, 1=user error (bad args/missing file), 2=data error (validation failed), 3=runtime error (OOM/failed computation)
- [x] **CLI-06**: Structured error reporting with error taxonomy (retryable, fixable, fatal) and actionable suggestion field
- [x] **CLI-07**: Non-interactive by default — no prompts, all params via flags/env/config. Fail fast on missing required params
- [x] **CLI-08**: Lazy imports — heavy dependencies (scanpy, torch, scvi-tools) loaded at command execution, not startup. `sct help` returns in <500ms

### CLI Core Commands

- [x] **CMD-01**: `sct qc run` — run QC metrics on a checkpoint file, output JSON metrics summary
- [x] **CMD-02**: `sct preprocess run` — run preprocessing recipe (modality-aware dispatch via `recipes.py`)
- [x] **CMD-03**: `sct validate <phase> <file>` — validate checkpoint against phase spec, return pass/fail JSON
- [x] **CMD-04**: `sct benchmark integration --from-dir <dir>` — load pre-computed h5ad files, compute comparison metrics, output ranked JSON
- [x] **CMD-05**: `sct status` — pipeline phase status from registry DAG (completed phases, available next, checkpoint paths)
- [x] **CMD-06**: `sct report <type>` — generate HTML report (pre-filter, post-filter, post-integration, post-celltyping)
- [x] **CMD-07**: Shared Result type used by both CLI commands and MCP tools (single implementation, dual serialization)
- [x] **CMD-08**: Fast-fail dependency check at command dispatch — report missing optional deps with install instructions before loading data

### CLI Discovery

- [ ] **DSC-01**: `sct list-commands --json` — machine-readable catalog of all commands with params, types, defaults
- [ ] **DSC-02**: `sct describe <cmd>` — JSON schema for specific command's params and output format
- [ ] **DSC-03**: `sct schema` — full CLI contract as JSON document (Typer command tree + Pydantic output schemas)

### Provenance

- [ ] **PRV-01**: JSON sidecar `.provenance.json` written alongside every CLI output file
- [ ] **PRV-02**: Sidecar includes: command, params, input files with SHA256 checksums, sc_tools version, timestamp, runtime_s, peak_memory_mb
- [ ] **PRV-03**: `sct provenance show <file>` — display provenance for a single output
- [ ] **PRV-04**: `sct provenance trace <file>` — trace full lineage DAG via input file references
- [ ] **PRV-05**: Reproducible Leiden clustering — configurable resolution and random_state propagation through all benchmark/clustering functions

### Scientific Gaps

- [ ] **SCI-01**: Pseudobulk DE module (`sc_tools.tl.de`) — aggregate counts by subject_id + celltype, run PyDESeq2 with batch covariates in design formula
- [ ] **SCI-02**: Marker validation report — after cell typing, compute top N marker genes per assigned type, generate dotplot, flag types with low canonical marker expression
- [ ] **SCI-03**: Subject-level metadata model — `subject_id` distinct from `library_id`, enforced at ingestion for multi-sample projects. Validate batch-condition confounding at registration
- [ ] **SCI-04**: Panel-aware cell typing dispatch — when `n_vars < 1000` (targeted panels), restrict to panel-validated methods and warn if whole-transcriptome model is applied

### Memory Safety

- [ ] **MEM-01**: IO Gateway with tiered loading strategy: h5py for metadata/embeddings, backed mode for summaries, full load only for compute
- [ ] **MEM-02**: `sct estimate <command> <args>` — pre-execution estimation of peak memory and runtime based on cell/gene count and method
- [ ] **MEM-03**: `--dry-run` flag for all data-touching commands — validate inputs, report what would happen, without executing

### Multi-Omic Assembly

- [ ] **MOM-01**: Patient/subject metadata join across modalities using validated `subject_id`
- [ ] **MOM-02**: MuData assembly from independently processed per-modality AnnData objects
- [ ] **MOM-03**: Joint embedding via multi-modal integration method (MOFA+, multiVI, or totalVI)
- [ ] **MOM-04**: Cross-modal queries at patient/project level (e.g., "show cell type proportions for patient X across all modalities")

### Testing

- [x] **TST-01**: Unit tests for `compute_integration_metrics` (synthetic data, single batch, no celltype)
- [x] **TST-02**: Unit tests for `compare_integrations` (NaN embeddings, single method, subsampling)
- [x] **TST-03**: Unit tests for `_stratified_subsample` (proportionality check, n > n_obs, single group)
- [x] **TST-04**: CLI argument parsing tests (no data loaded, fast)
- [x] **TST-05**: CLI integration tests with small fixtures (100-cell AnnData)
- [x] **TST-06**: End-to-end test with real data on HPC (skipif guard)

## v2 Requirements

Deferred to future milestone. Tracked but not in current roadmap.

### Advanced Analysis

- **ADV-01**: Spatial domain identification beyond UTAG (BankSY, GraphST)
- **ADV-02**: Ligand-receptor / cell-cell communication (CellChat, LIANA)
- **ADV-03**: Trajectory / RNA velocity analysis
- **ADV-04**: Reference atlas projection (scArches label transfer from CZ CELLxGENE Census)
- **ADV-05**: Deconvolution benchmarking (compare cell2location, tangram, destvi)
- **ADV-06**: Cell-typing benchmarking (compare sctype, celltypist, scArches systematically)
- **ADV-07**: Bootstrap uncertainty quantification for integration benchmarks (n_repeats parameter)
- **ADV-08**: Per-batch and per-celltype metric breakdown in integration benchmarks

### Infrastructure

- **INF-01**: W3C PROV-JSON export from sidecar provenance files
- **INF-02**: RO-Crate packaging for publication-ready reproducibility bundles
- **INF-03**: Database-backed provenance (migrate from sidecars once schema stabilizes)
- **INF-04**: CLI versioning / semver for output schema contracts

## Out of Scope

| Feature | Reason |
|---------|--------|
| Snakemake replacement | CLI complements workflow managers, doesn't replace them |
| Replacing MCP server | Both interfaces coexist, sharing the same backend |
| Web UI / dashboard | CLI-only; agents and scripts consume JSON |
| New integration methods | Use existing sc_tools implementations |
| Spatial coordinate registration | Project-specific; too variable for a CLI command |
| Daemon mode for import amortization | Assess need after Phase 2 startup benchmarks |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| BM-01 | Phase 1 | Complete |
| BM-02 | Phase 1 | Complete |
| BM-03 | Phase 1 | Complete |
| BM-04 | Phase 1 | Complete |
| BM-05 | Phase 1 | Complete |
| BM-06 | Phase 1 | Complete |
| BM-07 | Phase 1 | Complete |
| CLI-01 | Phase 2 | Complete |
| CLI-02 | Phase 2 | Complete |
| CLI-03 | Phase 2 | Complete |
| CLI-04 | Phase 2 | Complete |
| CLI-05 | Phase 2 | Complete |
| CLI-06 | Phase 2 | Complete |
| CLI-07 | Phase 2 | Complete |
| CLI-08 | Phase 2 | Complete |
| CMD-01 | Phase 3 | Complete |
| CMD-02 | Phase 3 | Complete |
| CMD-03 | Phase 3 | Complete |
| CMD-04 | Phase 3 | Complete |
| CMD-05 | Phase 3 | Complete |
| CMD-06 | Phase 3 | Complete |
| CMD-07 | Phase 3 | Complete |
| CMD-08 | Phase 3 | Complete |
| DSC-01 | Phase 4 | Pending |
| DSC-02 | Phase 4 | Pending |
| DSC-03 | Phase 4 | Pending |
| PRV-01 | Phase 5 | Pending |
| PRV-02 | Phase 5 | Pending |
| PRV-03 | Phase 5 | Pending |
| PRV-04 | Phase 5 | Pending |
| PRV-05 | Phase 5 | Pending |
| SCI-01 | Phase 6 | Pending |
| SCI-02 | Phase 6 | Pending |
| SCI-03 | Phase 6 | Pending |
| SCI-04 | Phase 6 | Pending |
| MEM-01 | Phase 7 | Pending |
| MEM-02 | Phase 7 | Pending |
| MEM-03 | Phase 7 | Pending |
| MOM-01 | Phase 8 | Pending |
| MOM-02 | Phase 8 | Pending |
| MOM-03 | Phase 8 | Pending |
| MOM-04 | Phase 8 | Pending |
| TST-01 | Phase 1 | Complete |
| TST-02 | Phase 1 | Complete |
| TST-03 | Phase 1 | Complete |
| TST-04 | Phase 2 | Complete |
| TST-05 | Phase 3 | Complete |
| TST-06 | Phase 3 | Complete |

**Coverage:**
- v1 requirements: 48 total
- Mapped to phases: 48
- Unmapped: 0 ✓

---
*Requirements defined: 2026-03-20*
*Last updated: 2026-03-20 after research + expert review*
