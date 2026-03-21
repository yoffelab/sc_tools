# Phase 2: CLI Foundation - Context

**Gathered:** 2026-03-21
**Status:** Ready for planning

<domain>
## Phase Boundary

`sct` is an installable, fast-starting CLI that produces structured JSON output and handles errors with actionable taxonomy. Requirements CLI-01 through CLI-08, TST-04. No actual pipeline commands yet — those come in Phase 3.

</domain>

<decisions>
## Implementation Decisions

### CLI module structure
- **D-01:** Start simple — single `sc_tools/cli.py` module, not a `sc_tools/cli/` package. Split into submodules later when Phase 3 adds command groups.
- **D-02:** `__main__.py` becomes a thin wrapper that imports and runs the Typer app. `python -m sc_tools` continues to work.
- **D-03:** Typer app object location is Claude's discretion (either in `cli.py` directly or a small `app.py` if import chains demand it).

### scripts/run_*.py deprecation
- **D-04:** Scripts (`run_qc_report.py`, `run_preprocessing.py`, `run_integration_benchmark.py`) are agent-written workarounds that bypass Snakemake. Phase 3 commands make them redundant, then they get deleted.
- **D-05:** No changes to scripts in Phase 2 — they stay as-is until `sct` commands replace them.

### CLIResult envelope
- **D-06:** Pydantic model with fields: `status`, `command`, `data`, `artifacts`, `provenance`, `message`.
- **D-07:** Status values: `success`, `error`, `partial`, `skipped`.
- **D-08:** `partial` status for operations that partially succeed (e.g., QC on 3/4 samples). Must include enough structure for an agent to identify and retry only the failed parts (parallelizable on cluster).
- **D-09:** `data` = computed results (metrics dict, status table). `artifacts` = file paths created (h5ad, HTML report, plots).
- **D-10:** Provenance field is minimal in Phase 2: `{command, timestamp, sc_tools_version}`. Full provenance (checksums, runtime, peak memory) deferred to Phase 5.
- **D-11:** Shared Result type introduced in Phase 2 — MCP tools also start returning it. Single implementation, dual serialization (JSON stdout for CLI, MCP return for tools).

### Error taxonomy & exit codes
- **D-12:** Top-level handler at the CLI boundary catches standard Python exceptions and maps them to the taxonomy. Additionally introduce sc_tools-specific exception classes (`SCToolsDataError`, `SCToolsRuntimeError`) that carry taxonomy metadata from the source.
- **D-13:** Three error categories: `retryable` (transient failures, OOM), `fixable` (user can correct input), `fatal` (unrecoverable).
- **D-14:** Actionable suggestions are templated per-error. Examples: fixable → `"Missing column 'batch'. Add batch annotation: adata.obs['batch'] = ..."`, retryable → `"OOM at 2.5M cells. Retry with --subsample-n 50000"`.
- **D-15:** Exit codes: 0=success (includes `partial`), 1=user error (bad args, missing file), 2=data error (file found but wrong format, validation failed), 3=runtime error (OOM, computation failure).
- **D-16:** `partial` results exit 0. Agent reads the JSON envelope to determine what succeeded/failed — exit codes are coarse pipeline signaling, the envelope carries the detail.

### Claude's Discretion
- Typer app object placement (cli.py vs app.py) based on import chain needs
- Lazy import implementation strategy (importlib, conditional imports, or Typer lazy groups)
- Rich output formatting specifics for `--human` mode
- Pydantic model details (field types, optional fields, serialization config)
- Test fixture design for CLI argument parsing tests (TST-04)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Existing CLI & MCP
- `sc_tools/__main__.py` — Current 40-line sys.argv parser; will become thin Typer wrapper
- `sc_tools/mcp/tools_server.py` — MCP server returning raw strings; will adopt shared Result type

### Package configuration
- `pyproject.toml` — No `[project.scripts]` entry yet; typer/rich/pydantic not in dependencies

### Pipeline phases (command groups will map to these)
- `sc_tools/pipeline.py` — `STANDARD_PHASES` with `PhaseSpec` objects; natural metadata for CLI command structure

### Requirements
- `.planning/REQUIREMENTS.md` — CLI-01 through CLI-08, TST-04 define exact acceptance criteria

### Prior phase context
- `.planning/phases/01-benchmark-fixes/01-CONTEXT.md` — Phase 1 decisions (no direct dependencies, but shared codebase patterns)

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `sc_tools/pipeline.py` `STANDARD_PHASES` — Phase names and specs can drive command group discovery
- `sc_tools/mcp/tools_server.py` — MCP tool patterns (validate_checkpoint, generate_qc_report) show what CLI commands will wrap
- `scripts/run_*.py` — Show the operations agents actually need (QC, preprocessing, benchmarking)

### Established Patterns
- MCP tools use `FastMCP` decorator pattern with typed parameters
- Heavy imports (scanpy, torch, scvi) are used throughout — lazy loading must be at the CLI layer, not deep in sc_tools
- `__main__.py` uses lazy `from sc_tools.registry import _cli_status` — precedent for deferred imports

### Integration Points
- `pyproject.toml` `[project.scripts]` — new `sct` entry point registration
- `pyproject.toml` dependencies — typer, rich, pydantic need to be added (pydantic may already be transitive)
- `__main__.py` — rewired to Typer app
- `sc_tools/mcp/tools_server.py` — adopt shared Result type for return values

</code_context>

<specifics>
## Specific Ideas

- `partial` results must be structured so agents can programmatically identify failed items and resubmit them as separate cluster jobs
- The scripts agents have been writing bypass Snakemake entirely — `sct` commands should be the replacement, not Snakemake integration
- Error suggestions should be specific enough that an agent can auto-fix without asking the user (e.g., exact flag to add, exact command to retry with)

</specifics>

<deferred>
## Deferred Ideas

- NaN guard in `pl/benchmarking.py` line 527 (`sc.pp.neighbors` on unfiltered embeddings) — not Phase 2 scope, note for future fix
- `scripts/run_*.py` deletion — Phase 3, after `sct` commands replace them
- Full provenance system — Phase 5
- Daemon mode for import amortization — assess after Phase 2 startup benchmarks

</deferred>

---

*Phase: 02-cli-foundation*
*Context gathered: 2026-03-21*
