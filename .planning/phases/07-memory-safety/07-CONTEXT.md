# Phase 7: Memory Safety - Context

**Gathered:** 2026-03-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Large datasets (2.5M cells, 25G h5ad) can be processed without OOM through tiered loading, pre-execution estimation, and dry-run validation. Requirements MEM-01, MEM-02, MEM-03. This phase adds an IO Gateway with tiered loading strategy, a `sct estimate` command for pre-execution memory/runtime prediction, and a `--dry-run` flag for all data-touching CLI commands.

</domain>

<decisions>
## Implementation Decisions

### IO Gateway location and scope
- **D-01:** New `sc_tools/io/` package — dedicated module with gateway logic. Separate from `storage.py` (which handles fsspec URI resolution). IO Gateway focuses on memory-aware loading strategy.
- **D-02:** CLI-path only — Gateway is used by `sct` CLI commands (via `cli_handler`). Library-level `sc_tools` functions keep calling `sc.read_h5ad` directly. Minimizes blast radius; agents get memory safety through CLI.

### Tier dispatch strategy
- **D-03:** Operation-based tier dispatch — each CLI command declares what access level it needs:
  - **T1 (h5py metadata-only):** Shape, column names, obsm keys. For `sct status`, `sct validate`, `sct estimate`, `sct provenance show`.
  - **T2 (backed/lazy mode):** obs/var DataFrames, obsm arrays without full X. For `sct qc run` (needs obs columns), `sct benchmark integration` (needs embeddings).
  - **T3 (full load):** Complete AnnData in memory. For `sct preprocess run`, `sct de run` (needs X matrix for computation).
- **D-04:** Gateway reads the command's tier declaration and loads accordingly. No size-based auto-promotion.

### Memory guard
- **D-05:** Pre-load memory guard with `--force` override. Gateway estimates memory from h5 metadata (n_obs × n_vars × dtype size) before T3 full load. If estimated memory exceeds available RAM, refuse with actionable error message suggesting `--force` flag or backed-mode alternative. Prevents the 44G OOM pattern that motivated this project.

### Claude's Discretion
- How `sct estimate` calculates memory/runtime projections (formula-based from cell×gene counts, method-specific profiles, or empirical lookup)
- Internal structure of tier declaration (decorator parameter, enum, or command metadata)
- How `--dry-run` is implemented (global flag on `cli_handler` intercepting before execution, or per-command logic)
- What `--dry-run` reports (input validation, planned operations, estimated resources, or all of the above)
- How backed mode (T2) interacts with operations that need partial X access (e.g., QC metrics on a subset of genes)
- Whether to reuse/extend `inspect_checkpoint()` h5py pattern from MCP or write fresh gateway code
- Memory estimation formula calibration (static multipliers vs empirical profiling)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Existing memory infrastructure
- `sc_tools/memory/profiling.py` — `get_memory_usage()`, `estimate_adata_memory()`, `check_memory_threshold()`, `aggressive_cleanup()`. Current estimation only works on loaded AnnData objects.
- `sc_tools/memory/__init__.py` — Module exports

### Existing IO patterns (reuse candidates)
- `sc_tools/storage.py` — `smart_read_h5ad()` with backed mode support, fsspec URI resolution
- `sc_tools/utils/checkpoint.py` — `read_checkpoint()` with `backed` parameter
- `sc_tools/bm/integration.py` §63-110 — `_load_embedding_h5py()` — proven h5py selective read pattern for embeddings
- `sc_tools/mcp/tools_server.py` §487-530 — `inspect_checkpoint()` — h5py metadata-only reader (shape, column names, keys)

### CLI infrastructure (integration points)
- `sc_tools/cli/__init__.py` §223+ — `cli_handler` decorator — interception point for --dry-run and memory guard
- `sc_tools/models/result.py` — `CLIResult` envelope with `peak_memory_mb` in provenance
- `sc_tools/provenance/sidecar.py` §32 — `get_peak_memory_mb()` via `resource.getrusage`

### Existing dry-run patterns
- `sc_tools/bm/slurm.py` §151-170 — `submit_benchmark_jobs()` dry_run parameter
- `sc_tools/bm/cli.py` §86 — `--dry-run` CLI argument for SLURM submission
- `sc_tools/mcp/tools_server.py` §706-824 — `run_full_phase()` dry_run parameter

### Prior phase decisions
- `sc_tools/provenance/sidecar.py` — Phase 5: peak_memory_mb already tracked in provenance sidecars
- `.planning/phases/05-provenance/05-CONTEXT.md` — Provenance model decisions (D-06 through D-08)

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `memory/profiling.py`: `estimate_adata_memory()` computes memory from X matrix + obs + var + obsm. Can be adapted for pre-load estimation from h5 metadata.
- `bm/integration.py`: `_load_embedding_h5py()` reads obsm arrays via h5py without full AnnData load — proven T1/T2 pattern.
- `mcp/tools_server.py`: `inspect_checkpoint()` reads h5ad metadata (shape, column names, keys) via h5py — proven T1 pattern.
- `storage.py`: `smart_read_h5ad()` already supports `backed` parameter — T2 loading exists.
- `provenance/sidecar.py`: `get_peak_memory_mb()` measures actual peak memory — useful for estimate calibration.

### Established Patterns
- `cli_handler` decorator is the standard interception point for all CLI commands — natural place for --dry-run and memory guard.
- Commands return `CLIResult` with data/artifacts/provenance — dry-run can return a CLIResult with status="dry_run" and estimated resources in data.
- Lazy imports via `_check_deps()` — IO Gateway should follow same pattern for h5py dependency.

### Integration Points
- `cli_handler` in `cli/__init__.py` — add --dry-run interception and memory guard check before command execution
- Each CLI command file (qc.py, preprocess.py, benchmark.py, de.py) — add tier declaration metadata
- `sct estimate` would be a new top-level command (like `sct provenance`) registered in `cli/__init__.py`

</code_context>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 07-memory-safety*
*Context gathered: 2026-03-24*
