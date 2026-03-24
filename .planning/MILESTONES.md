# Milestones

## v1.0 Agent-Native CLI (Shipped: 2026-03-24)

**Phases completed:** 8 phases, 18 plans, 33 tasks

**Key accomplishments:**

- Fixed five integration benchmark bugs (NaN handling, subsampling bias, subsample_n, runtime tracking, param provenance) with 15 new TDD tests covering edge cases
- h5py-based embedding loading for 2.5M-cell benchmarks without OOM, bio_key parameter for flexible bio conservation, and fixed _recipe_targeted_panel to preserve raw counts for scVI
- CLIResult Pydantic envelope with Status/ErrorInfo/Provenance models and 5-class exception hierarchy mapping to semantic exit codes 0/1/2/3
- Typer CLI app with global --human flag, error handler mapping exceptions to exit codes 0-3, five stub command groups, 24 argument parsing tests, and MCP CLIResult migration proof-of-concept
- cli.py migrated to cli/ package with _check_deps fast-fail utility and adata_100 shared test fixture
- Four CLI commands wrapping existing backend: validate (checkpoint checks), status (DAG viewer with registry fallback), qc run (per-sample metrics), and report generate (HTML report dispatch)
- sct preprocess run with modality auto-detection (D-08) and sct benchmark integration with --from-dir embedding discovery and --report HTML generation (D-12)
- 23 integration tests covering all CLI commands with CLIResult JSON verification, E2E scaffold with skipif guards, and verified CLI/MCP shared Result type
- Three discovery commands (list-commands, describe, schema) providing machine-readable CLI introspection via Typer/Click tree walking and Pydantic model_json_schema
- ProvenanceRecord model with SHA256 checksums, automatic .provenance.json sidecar writing in cli_handler, h5ad uns embedding via h5py, and random_state threading through Leiden clustering
- BFS lineage trace engine with cycle detection, SHA256 relocation, and adata.uns fallback; sct provenance show/trace CLI commands registered via register_provenance(app) pattern
- Subject-level metadata validation with confounding detection, and panel-aware cell typing dispatch restricting whole-transcriptome methods for targeted panels (n_vars < 1000)
- Pseudobulk DE module wrapping PyDESeq2 with per-celltype aggregation, collinearity-guarded design formula, and sct de run CLI command
- Marker validation compute + post-celltyping report integration with threshold-based flagging, dotplot, and flag table
- Tiered IO Gateway with h5py metadata reader, memory estimator, and 80%-RAM guard for safe large-file loading
- sct estimate command with memory/runtime projection and --dry-run/--force flags on all data-touching CLI commands
- MultiOmicAtlas class with outer-join patient metadata, MuData construction, and cross-modal celltype proportion queries
- 1. [Rule 1 - Bug] Return type annotation causes Typer NameError on Python 3.10

---
