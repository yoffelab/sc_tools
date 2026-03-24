# Phase 3: Core Commands - Context

**Gathered:** 2026-03-21
**Status:** Ready for planning

<domain>
## Phase Boundary

Agents can run QC, preprocessing, validation, benchmarking, status checks, and report generation through `sct` commands instead of ad-hoc scripts. Requirements CMD-01 through CMD-08, TST-05, TST-06. This phase fills the 5 stub command groups created in Phase 2 with real implementations. After this phase, `scripts/run_*.py` become redundant and can be deleted.

</domain>

<decisions>
## Implementation Decisions

### CLI module structure
- **D-01:** Split `cli.py` into a `cli/` package. Current `cli.py` content (app, cli_handler, shared utils) moves to `cli/__init__.py`. Commands go into per-domain submodules:
  - `cli/qc.py` — `sct qc run`, `sct report`
  - `cli/preprocess.py` — `sct preprocess run`
  - `cli/benchmark.py` — `sct benchmark integration`
  - `cli/validate.py` — `sct validate`
  - `cli/status.py` — `sct status`
- **D-02:** No external API change — `from sc_tools.cli import app` continues to work.

### Output path conventions
- **D-03:** All commands support `--project-dir` flag (default `.`). Output paths are relative to project dir.
- **D-04:** Report output: `{project_dir}/figures/reports/{YYMMDD}_{input_file}_{report_type}_report.html`. Figures: `{project_dir}/figures/reports/{YYMMDD}_figures/`. No `--figures-dir` flag needed — auto-derived from report path.
- **D-05:** Report type is always HTML. No `--report-type` flag.

### Command argument design (sensible defaults, all overridable)
- **D-06:** Every flag has a sensible default. Commands should "just work" with minimal args for the common case. All parameters remain available as flags for override.
- **D-07:** `sct benchmark integration --from-dir <dir>`:
  - `--subsample-n 500000` (default; auto-applied when dataset > 500k cells)
  - `--seed 0` (default)
  - When subsampling occurs, log clearly: "Subsampling from {N} to 500k cells (seed=0)"
  - Record subsampling decision in CLIResult provenance AND in benchmark report HTML
  - Benchmark report must clearly communicate per-sample how results were generated to avoid confusion
- **D-08:** `sct preprocess run`:
  - `--modality auto` (default) — auto-detect from `adata.uns` metadata
  - Throw `SCToolsDataError` if modality cannot be determined (with suggestion: "Specify --modality visium|visium_hd|xenium|cosmx|imc")
- **D-09:** `sct validate <phase> <file>` and `sct status`:
  - Graceful fallback when registry is unavailable — return what info is available without erroring
  - `sct status` without registry shows phase DAG from `pipeline.py` with "registry unavailable" note

### Shared Result type (CMD-07)
- **D-10:** All CLI commands return CLIResult. MCP tools that wrap the same operations share the backend function and return CLIResult. Phase 2 PoC (validate_checkpoint) already demonstrates this pattern — extend to all commands.
- **D-11:** Fast-fail dependency check (CMD-08): before loading data, check if required optional deps are available. Report missing deps with install instructions (e.g., "scvi-tools required for scVI integration: pip install scvi-tools").

### Transparency & communication
- **D-12:** When benchmark subsamples, the decision must appear in:
  1. CLI log output (always)
  2. CLIResult provenance field
  3. Benchmark HTML report (per-sample section)
- **D-13:** All commands that modify or generate data should communicate clearly what was done, what parameters were used, and how to reproduce.

### scripts/run_*.py deprecation
- **D-14:** After Phase 3 commands are verified working, delete `scripts/run_qc_report.py`, `scripts/run_preprocessing.py`, `scripts/run_integration_benchmark.py`. Keep `scripts/run_qc_reports_all_projects.py` if it has unique multi-project logic not covered by the CLI.

### Claude's Discretion
- Exact flag names and short aliases for each command
- Output file naming details beyond the `{YYMMDD}_{input}_{type}` pattern
- How `sct report` dispatches to different report types (pre_filter, post_filter, post_integration, post_celltyping)
- `sct status` output format (table vs tree vs list)
- Internal organization of CLI submodules (helper functions, shared utilities)
- How `sct qc run` structures its JSON metrics summary

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### CLI foundation (Phase 2 output)
- `sc_tools/cli.py` — Current Typer app with 5 stub groups, cli_handler, --human flag (will become cli/__init__.py)
- `sc_tools/models/result.py` — CLIResult, Status, ErrorInfo, Provenance models
- `sc_tools/errors.py` — Exception hierarchy (SCToolsUserError, SCToolsDataError, SCToolsRuntimeError, SCToolsFatalError)

### Backend functions to wrap
- `sc_tools/qc/sample_qc.py` — QC metrics computation
- `sc_tools/qc/report.py` — HTML report generation (pre_filter, post_filter, post_integration)
- `sc_tools/pp/recipes.py` — Modality-aware preprocessing dispatch (`_recipe_targeted_panel`, `_recipe_standard`)
- `sc_tools/bm/integration.py` — `compare_integrations()`, `compute_integration_metrics()`, `run_integration_benchmark()`
- `sc_tools/validate.py` — Checkpoint validation against PhaseSpec
- `sc_tools/pipeline.py` — `STANDARD_PHASES`, `get_available_next()`, `get_phase_checkpoint()`

### MCP tools (dual serialization targets)
- `sc_tools/mcp/tools_server.py` — `validate_checkpoint` (already returns CLIResult), `generate_qc_report`, `load_sample`, `run_full_phase`

### Existing scripts (reference for argument patterns)
- `scripts/run_qc_report.py` — QC report workflow with --report, --adata, --adata-post flags
- `scripts/run_preprocessing.py` — Preprocessing script
- `scripts/run_integration_benchmark.py` — Integration benchmark script

### Requirements
- `.planning/REQUIREMENTS.md` — CMD-01 through CMD-08, TST-05, TST-06

### Prior context
- `.planning/phases/02-cli-foundation/02-CONTEXT.md` — CLIResult design, error taxonomy, all Phase 2 decisions

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `cli_handler` decorator in `cli.py` — wraps any command function, handles errors, emits CLIResult
- `CLIResult.model_dump_json()` for CLI, `.model_dump(mode="json")` for MCP — dual serialization
- `_emit()` and `_render_rich()` in `cli.py` — stdout/stderr output helpers
- `compare_integrations()` already supports `subsample_n`, `seed`, `embedding_files` (Phase 1 work)
- `recipes.py` dispatches by modality from `adata.uns`

### Established Patterns
- Lazy imports inside command functions (Phase 2, CLI-08)
- `@cli_handler` decorator on every command function
- Typer `add_typer()` for command groups with stub callbacks

### Integration Points
- `cli.py` → `cli/__init__.py` migration (move app + shared utils)
- Each `cli/*.py` submodule registers commands on its group's Typer app
- MCP tools call the same backend functions as CLI commands
- `conftest.py` in `sc_tools/tests/` needs creation for shared fixtures

</code_context>

<specifics>
## Specific Ideas

- Benchmark report HTML should clearly show per-sample: original cell count, whether subsampled, subsample_n, seed, methods compared — so a user coming back months later understands how results were generated
- Auto-detection for modality (`--modality auto`) should check `adata.uns` first, fall back to checking panel size (`n_vars < 1000` suggests targeted panel)
- `sct status` should be useful even without a registry — show the phase DAG from `pipeline.py` and note which phases have checkpoint files on disk

</specifics>

<deferred>
## Deferred Ideas

- NaN guard in `pl/benchmarking.py` line 527 — carried forward from Phase 2
- Daemon mode for import amortization — assess after Phase 3
- `sct celltype` commands — stub group exists but no commands in Phase 3 requirements
- Multi-project batch commands (like `run_qc_reports_all_projects.py`) — future phase if needed

</deferred>

---

*Phase: 03-core-commands*
*Context gathered: 2026-03-21*
