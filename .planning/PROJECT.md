# sc_tools Agent-Native CLI

## What This Is

A CLI-first interface layer on top of sc_tools that makes computational biology pipelines agent-native. Instead of agents writing one-off scripts when sc_tools has interface gaps, they call stable, composable CLI commands with structured output. The tool covers the full single-cell/spatial analysis lifecycle: ingest, QC, preprocessing, integration, cell typing, benchmarking, and multi-omic assembly — with lightweight provenance tracking throughout.

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
- ✓ Cell type annotation (sctype, celltypist, scArches) — existing `sc_tools/tl/`
- ✓ Integration benchmarking metrics (silhouette, cLISI, iLISI, kBET) — existing `sc_tools/bm/`
- ✓ MCP server exposing core operations — existing `sc_tools/mcp/`
- ✓ GPU-aware execution with CPU fallback — existing `sc_tools/pp/_gpu.py`
- ✓ Registry database for phase/provenance tracking — existing `sc_tools/migrations/`

### Active

<!-- Current scope. Building toward these. -->

- [ ] Pre-computed integration benchmark support (load separate per-method h5ad files via h5py, NaN handling, subsampling)
- [ ] Agent-native CLI interface wrapping sc_tools operations (verb-based commands, JSON stdout, self-describing help)
- [ ] Structured output format for all CLI commands (machine-readable JSON that agents parse to decide next steps)
- [ ] CLI self-discovery (`sct help`, `sct list-commands`) so agents know what's available
- [ ] Lightweight file-based provenance tracking (JSON sidecar files alongside outputs)
- [ ] Patient/subject metadata attachment post-ingestion
- [ ] MuData assembly from independently processed modalities
- [ ] Cross-modal multi-omic queries at patient/project level

### Out of Scope

- Full Snakemake replacement — CLI complements Snakemake, doesn't replace it
- Web UI / dashboard — CLI and structured output only
- New integration methods — use existing sc_tools implementations
- Database-first provenance schema — start file-based, migrate to DB once model stabilizes
- Changing run_preprocessing.py to save only embeddings — separate optimization

## Context

**Origin:** Agents (Claude Code) running comp bio workflows default to writing one-off scripts when sc_tools has interface gaps. Example: integration benchmarking ran Harmony/scVI/resolVI in parallel (good), but sc_tools couldn't consume the separate outputs into a benchmark report. The agent wrote a custom script that OOM'd at 44G. This pattern repeats across preprocessing, filtering, and cell typing.

**Existing codebase:** sc_tools is a mature library (Python 3.11+, scanpy-based, AnnData-centric) with ingestion, QC, preprocessing, integration, analysis, and benchmarking modules. It has an MCP server and a registry database, but the CLI interface is ad-hoc (scattered scripts in `scripts/`).

**Data scale:** Visium HD datasets reach 2.5M cells, 11-25G per h5ad file. Memory-efficient access (h5py, backed mode) is essential.

**Registry status:** SQLAlchemy-based registry exists but schema is unstable — frequent migrations, weak typing. The schema needs to emerge from actual usage before committing to a formal data model.

**Multi-modal landscape:** scRNA-seq, IMC, Visium, Visium HD, Xenium — different raw formats, different processing pipelines, all eventually tied to patients/projects for multi-omic questions.

## Constraints

- **Engine:** sc_tools is the backend — CLI wraps it, doesn't replace it
- **Memory:** Must handle 2.5M-cell datasets without OOM — h5py for selective reads, subsampling for metrics
- **Agent compatibility:** CLI must work for both agents (JSON output) and humans (readable output with `--human` or similar flag)
- **Provenance:** File-based first (JSON sidecars), DB later once domain model stabilizes
- **Platforms:** Linux (HPC: BRB, Cayuga) and macOS (local dev)

## Key Decisions

<!-- Decisions that constrain future work. Add throughout project lifecycle. -->

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| CLI wraps sc_tools, not replaces it | Existing internals work; the gap is the interface layer | — Pending |
| File-based provenance before DB | Schema instability means DB migrations will keep breaking; let model emerge from usage | — Pending |
| Fix sc_tools gaps first, then build CLI | Immediate value (integration benchmark), and CLI needs stable internals to wrap | — Pending |
| MuData assembly is late-stage | Each modality processed independently; assembly happens after analysis, not at ingestion | — Pending |
| JSON structured output by default | Agent-native design; human-readable mode as opt-in flag | — Pending |

---
*Last updated: 2026-03-20 after initialization*
