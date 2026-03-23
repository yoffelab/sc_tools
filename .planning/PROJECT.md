# sc_tools Agent-Native CLI

## What This Is

A CLI-first interface layer (`sct`) on top of sc_tools that makes computational biology pipelines agent-native. Instead of agents writing one-off scripts when sc_tools has interface gaps, they call stable, composable CLI commands with structured JSON output. Built on Typer with a Pydantic CLIResult envelope, the tool covers the full single-cell/spatial analysis lifecycle — ingest, QC, preprocessing, integration, benchmarking, cell typing, differential expression, validation, and multi-omic assembly — with file-based provenance tracking throughout.

The CLI coexists with the existing MCP server. Both interfaces share the same backend functions via a common Result type. CLI is for subprocess invocation by agents; MCP is for in-process tool calling.

## Core Value

Agents never write throwaway scripts. Every comp bio operation is callable via a stable CLI with structured I/O, so analysis is reproducible by default.

## Requirements

### Validated

<!-- Shipped and confirmed valuable — existing sc_tools capabilities. -->

- ✓ Platform-specific data ingestion (Visium, Visium HD, Xenium, CosMx, IMC) — existing `sc_tools/ingest/`
- ✓ QC metrics calculation and sample classification — existing `sc_tools/qc/`
- ✓ HTML report generation (pre-filter, post-filter, post-integration, post-celltyping) — existing `sc_tools/qc/report.py`
- ✓ Normalization, HVG selection, PCA, UMAP, clustering — existing `sc_tools/pp/`
- ✓ Batch integration (scVI, Harmony, BBKNN, scanorama) — existing `sc_tools/pp/`
- ✓ Gene set scoring and enrichment analysis — existing `sc_tools/tl/`
- ✓ Cell type annotation (sctype, celltypist, ensemble) — existing `sc_tools/tl/`
- ✓ Integration benchmarking metrics (silhouette, cLISI, iLISI, kBET) — existing `sc_tools/bm/`
- ✓ MCP server exposing core operations — existing `sc_tools/mcp/`
- ✓ GPU-aware execution with CPU fallback — existing `sc_tools/pp/_gpu.py`
- ✓ Registry database for phase/provenance tracking — existing `sc_tools/migrations/`
- ✓ Pipeline phase DAG with `PhaseSpec` metadata — existing `sc_tools/pipeline.py`
- ✓ Modality-aware preprocessing recipes — existing `sc_tools/pp/recipes.py`

### Active

<!-- Current scope. Building toward these. -->

**sc_tools internals fixes (Phase 1) — Validated 2026-03-21:**
- [x] Pre-computed integration benchmark loading (h5py selective reads, NaN filtering per-embedding, configurable subsampling)
- [x] Fix `_stratified_subsample` truncation bias (proportional downsampling, not `sorted()[:n]`)
- [x] Benchmark parameter provenance (batch_weight, bio_weight, seed, resolution stored in output)
- [x] Runtime tracking in integration benchmarks (`runtime_s` column)
- [x] NaN handling in `compare_integrations()` (per-embedding valid-cell masking)
- [x] Fix `_recipe_targeted_panel` normalizing before scVI raw count check

**CLI foundation (Phase 2) — Validated 2026-03-21:**
- [x] `sct` entry point registered in `pyproject.toml` `[project.scripts]`
- [x] Typer-based CLI framework with Pydantic CLIResult envelope (`{status, command, data, artifacts, provenance, message}`)
- [x] JSON stdout by default, `--human` flag for Rich-formatted stderr output
- [x] Semantic exit codes (0=success, 1=user error, 2=data error, 3=runtime error)
- [x] Structured error reporting with retryable/fixable/fatal taxonomy
- [x] Non-interactive by default (no prompts; all params via flags/env/config)
- [x] Lazy imports (heavy deps like scanpy/torch loaded at command execution, not startup; <500ms for `sct help`)

**CLI core commands (Phase 3) — Validated 2026-03-22:**
- [x] `sct qc run`, `sct preprocess run`, `sct validate <phase>` — wrap existing operations
- [x] `sct benchmark integration --from-dir <dir>` — pre-computed benchmark comparison
- [x] `sct status` — pipeline phase status from registry/DAG
- [x] Shared Result type used by both CLI and MCP tools
- [x] Fast-fail for missing optional dependencies at command dispatch

**CLI self-discovery (Phase 4):**
- [ ] `sct list-commands --json` — machine-readable command catalog
- [ ] `sct describe <cmd>` — JSON schema for command params and output
- [ ] `sct schema` — full CLI contract as JSON (Typer introspection + Pydantic schemas)

**Provenance & reproducibility (Phase 5):**
- [ ] JSON sidecar provenance files (`.provenance.json` alongside every output)
- [ ] Provenance includes: command, params, input files with checksums, sc_tools version, timestamp, runtime, peak memory
- [ ] `sct provenance show <file>` and `sct provenance trace <file>` for lineage queries
- [ ] Reproducible Leiden resolution (configurable, with random_state propagation)

**Scientific gaps (Phase 6):**
- [ ] Pseudobulk DE module (`sc_tools.tl.de`) wrapping PyDESeq2 with automatic aggregation by subject_id + celltype
- [ ] Marker validation report after cell typing (dotplot/heatmap of canonical markers per assigned type)
- [ ] Subject-level metadata model (`subject_id` distinct from `library_id`, enforced at ingestion for multi-sample projects)
- [ ] Panel-aware cell typing dispatch (restrict methods when `n_vars < 1000`)

**Memory safety (Phase 7):**
- [ ] IO Gateway with tiered loading strategy (h5py for metadata, backed for summaries, full for compute)
- [ ] `sct estimate <command>` — pre-execution memory/time estimation
- [ ] Dry-run mode for all data-touching commands

**Multi-omic assembly (Phase 8):**
- [ ] Patient/subject metadata join across modalities
- [ ] MuData assembly from independently processed per-modality AnnData
- [ ] Joint embedding via multi-modal integration (MOFA+, multiVI)
- [ ] Cross-modal queries at patient/project level

### Out of Scope

- Full Snakemake replacement — CLI complements Snakemake, doesn't replace it
- Replacing MCP server with CLI — both interfaces coexist, sharing the same backend
- Web UI / dashboard — CLI and structured output only
- New integration methods — use existing sc_tools implementations
- Database-first provenance schema — start file-based, migrate to DB once model stabilizes
- Spatial coordinate registration across platforms — project-specific, not CLI command
- Deconvolution benchmarking — use external tools until demand proven
- Cell-typing benchmarking — defer until annotation methods are more stable
- Trajectory / RNA velocity — not in current pipeline, defer to future milestone
- Ligand-receptor / cell-cell communication — defer to future milestone
- Spatial domain identification beyond UTAG — defer to future milestone
- Daemon mode for import amortization — assess if needed after Phase 2

## Context

**Origin:** Agents (Claude Code) running comp bio workflows default to writing one-off scripts when sc_tools has interface gaps. Example: integration benchmarking ran Harmony/scVI/resolVI in parallel (good), but sc_tools couldn't consume the separate outputs into a benchmark report. The agent wrote a custom script that OOM'd at 44G. This pattern repeats across preprocessing, filtering, and cell typing.

**Existing codebase:** sc_tools is a mature library (Python 3.11+, scanpy-based, AnnData-centric) with ingestion, QC, preprocessing, integration, analysis, and benchmarking modules. It has an MCP server (8+ tools) and a registry database, but the CLI interface is a 40-line `sys.argv` parser in `__main__.py` handling one command (`registry status`). No `sct` entry point is registered in `pyproject.toml`. The `scripts/run_*.py` collection is the de facto CLI.

**Pipeline DAG:** `sc_tools/pipeline.py` defines `STANDARD_PHASES` with `PhaseSpec` objects including checkpoint paths, required obs/obsm, dependencies, and human-in-loop flags. This is natural metadata for CLI self-description.

**Data scale:** Visium HD datasets reach 2.5M cells, 11-25G per h5ad file. Full load expands to ~4x in memory (~100G). Memory-efficient access (h5py, backed mode) is essential. The 44G OOM that motivated this project is a recurring pattern.

**Registry status:** SQLAlchemy-based registry with 18 migrations, 7 unapplied (0012-0018). Four-layer schema (Projects → Datasets → BioData → DataSources). Not unstable in code — unstable in deployment. File-based provenance sidecars avoid the migration treadmill while the model stabilizes.

**Multi-modal landscape:** scRNA-seq, IMC, Visium, Visium HD, Xenium — different raw formats, different processing pipelines, all eventually tied to patients/projects for multi-omic questions.

**Known code bugs (from review):**
- `_stratified_subsample`: `sorted(indices)[:n]` truncation biases toward low array indices
- `_recipe_targeted_panel`: normalizes data before scVI raw count check (data-destructive)
- `extract_reference_profiles`: negative shift distorts relative expression differences
- `_asw_batch_sklearn` vs `scib_metrics.silhouette_batch`: same column name, different computations, no indicator

## Constraints

- **Engine:** sc_tools is the backend — CLI wraps it, doesn't replace it
- **Memory:** Must handle 2.5M-cell datasets without OOM — tiered loading strategy (h5py/backed/full)
- **Agent compatibility:** JSON to stdout, logging/progress to stderr — never mix. `--human` flag for Rich output
- **Startup performance:** `sct help` must return in <500ms. No heavy imports (scanpy, torch) at startup
- **Provenance:** File-based first (JSON sidecars), DB later once domain model stabilizes
- **CLI stability:** Command names and JSON output shapes are API contracts. Breaking changes need versioning
- **Benchmark reproducibility:** All metrics that depend on clustering must document resolution and random state
- **Pseudobulk DE over single-cell DE:** For multi-sample comparisons, pseudobulk with proper sample-level modeling (Squair et al. 2021)
- **Platforms:** Linux (HPC: BRB, Cayuga) and macOS (local dev)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| CLI wraps sc_tools, not replaces it | Existing internals work; the gap is the interface layer | — Pending |
| CLI uses Typer for argument parsing | Type-hint-driven help matches codebase conventions; auto-completion support | — Pending |
| CLI and MCP share a Result return type | Avoid reimplementing structured output twice; single source of truth | — Pending |
| CLI commands map to pipeline phases, not Python modules | Agents think in pipeline steps (`sct qc`), not library imports (`sct qc.metrics`) | — Pending |
| Register `sct` in pyproject.toml `[project.scripts]` | Makes CLI installable and discoverable via `which sct` | — Pending |
| JSON structured output by default | Agent-native design; human-readable mode as opt-in flag via `--human` | — Pending |
| File-based provenance before DB | Schema instability means DB migrations will keep breaking; let model emerge from usage | — Pending |
| Fix sc_tools gaps first, then build CLI | Immediate value (integration benchmark), and CLI needs stable internals to wrap | — Pending |
| MuData assembly is late-stage | Each modality processed independently; assembly happens after analysis, not at ingestion | — Pending |
| Late-stage MuData uses patient-level join, not spatial registration | Spatial co-registration is project-specific; CLI supports metadata joins robustly | — Pending |
| Pseudobulk DE over single-cell DE for multi-sample comparisons | Single-cell tests inflate significance due to pseudoreplication | — Pending |

---
*Last updated: 2026-03-22 — Phase 3 (core commands) complete*
