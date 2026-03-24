# Phase 3: Core Commands - Research

**Researched:** 2026-03-21
**Domain:** Typer CLI command implementation wrapping existing sc_tools backend functions
**Confidence:** HIGH

## Summary

Phase 3 fills the five stub command groups created in Phase 2 with real implementations. The work is primarily a wiring exercise: each CLI command wraps an existing backend function (QC, preprocessing, validation, benchmarking, status) with Typer argument definitions, the `@cli_handler` decorator, and CLIResult envelope construction. The `cli.py` monolithic module must be split into a `cli/` package per decision D-01.

The existing backend functions (`sample_qc.py`, `recipes.py`, `integration.py`, `validate.py`, `pipeline.py`, `report.py`) are mature and well-tested. The existing scripts (`run_qc_report.py`, `run_preprocessing.py`, `run_integration_benchmark.py`) provide the exact argument patterns and data-loading boilerplate that each CLI command must replicate. The MCP server (`tools_server.py`) already demonstrates the CLIResult dual-serialization pattern with `validate_checkpoint`.

**Primary recommendation:** Structure work as cli/ package migration first (Wave 0), then implement commands incrementally by complexity (validate/status first, then qc/report, then preprocess, then benchmark), with CMD-07/CMD-08 as cross-cutting concerns woven into every command.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Split `cli.py` into a `cli/` package. Current `cli.py` content (app, cli_handler, shared utils) moves to `cli/__init__.py`. Commands go into per-domain submodules: `cli/qc.py`, `cli/preprocess.py`, `cli/benchmark.py`, `cli/validate.py`, `cli/status.py`
- **D-02:** No external API change -- `from sc_tools.cli import app` continues to work
- **D-03:** All commands support `--project-dir` flag (default `.`). Output paths are relative to project dir
- **D-04:** Report output: `{project_dir}/figures/reports/{YYMMDD}_{input_file}_{report_type}_report.html`
- **D-05:** Report type is always HTML. No `--report-type` flag
- **D-06:** Every flag has a sensible default. Commands should "just work" with minimal args
- **D-07:** `sct benchmark integration --from-dir <dir>`: `--subsample-n 500000` default, `--seed 0` default, with logging and provenance recording of subsampling decisions
- **D-08:** `sct preprocess run`: `--modality auto` default, throw `SCToolsDataError` if modality cannot be determined
- **D-09:** `sct validate` and `sct status`: graceful fallback when registry unavailable
- **D-10:** All CLI commands return CLIResult. MCP tools share the same backend function
- **D-11:** Fast-fail dependency check before loading data
- **D-12:** Subsampling decisions appear in CLI log, CLIResult provenance, and benchmark HTML report
- **D-13:** All commands communicate clearly what was done, parameters used, and how to reproduce
- **D-14:** After Phase 3, delete `scripts/run_qc_report.py`, `scripts/run_preprocessing.py`, `scripts/run_integration_benchmark.py`

### Claude's Discretion
- Exact flag names and short aliases for each command
- Output file naming details beyond the `{YYMMDD}_{input}_{type}` pattern
- How `sct report` dispatches to different report types (pre_filter, post_filter, post_integration, post_celltyping)
- `sct status` output format (table vs tree vs list)
- Internal organization of CLI submodules (helper functions, shared utilities)
- How `sct qc run` structures its JSON metrics summary

### Deferred Ideas (OUT OF SCOPE)
- NaN guard in `pl/benchmarking.py` line 527 -- carried forward from Phase 2
- Daemon mode for import amortization -- assess after Phase 3
- `sct celltype` commands -- stub group exists but no commands in Phase 3 requirements
- Multi-project batch commands -- future phase if needed
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| CMD-01 | `sct qc run` -- run QC metrics on checkpoint, output JSON summary | Wraps `compute_sample_metrics()` + `classify_samples()` from `sample_qc.py`. Pattern from `run_qc_report.py` lines 109-122 |
| CMD-02 | `sct preprocess run` -- modality-aware preprocessing | Wraps `preprocess()` from `recipes.py`. Auto-detect modality from `adata.uns` or panel size |
| CMD-03 | `sct validate <phase> <file>` -- checkpoint validation | Wraps `validate_file()` from `validate.py`. MCP pattern already in `tools_server.py` |
| CMD-04 | `sct benchmark integration --from-dir <dir>` -- benchmark comparison | Wraps `compare_integrations()` from `integration.py` with `embedding_files` mode |
| CMD-05 | `sct status` -- pipeline DAG state | Uses `get_dag()`, `get_available_next()`, `get_phase_checkpoint()` from `pipeline.py` |
| CMD-06 | `sct report <type>` -- HTML report generation | Dispatches to `generate_{pre_filter,post_filter,post_integration,post_celltyping}_report()` |
| CMD-07 | Shared Result type for CLI and MCP | CLIResult already exists; extend MCP tools to call same backend functions |
| CMD-08 | Fast-fail dependency check | Check for optional deps (scanpy, scvi-tools, etc.) before loading data |
| TST-05 | CLI integration tests with 100-cell AnnData fixtures | Need `conftest.py` with shared fixtures, CliRunner-based tests |
| TST-06 | End-to-end test with real data on HPC | `skipif` guard when HPC/real data not available |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| typer | 0.24.1 | CLI framework | Already in use (Phase 2), sub-app pattern established |
| pydantic | 2.x | CLIResult model | Already in use (Phase 2), dual serialization proven |
| rich | installed | Human-readable output | Already wired via `_render_rich()` in cli.py |
| anndata | installed | AnnData I/O | Core data format for all commands |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| h5py | installed | Lightweight embedding reads | `sct benchmark` reads obsm without full AnnData load |
| scanpy | installed | QC metrics, preprocessing | Lazy-imported inside command functions |
| pytest + typer.testing.CliRunner | installed | CLI integration tests | TST-05 test infrastructure |

### Alternatives Considered
None -- all libraries are locked by Phase 2 decisions. No new dependencies needed.

## Architecture Patterns

### Recommended Project Structure
```
sc_tools/
  cli/
    __init__.py     # app, _state, main(), _emit(), _render_rich(), cli_handler,
                    # _make_error_result, _check_deps(), sub-app registrations
    qc.py           # sct qc run (CMD-01), sct report (CMD-06)
    preprocess.py   # sct preprocess run (CMD-02)
    validate.py     # sct validate <phase> <file> (CMD-03)
    benchmark.py    # sct benchmark integration (CMD-04)
    status.py       # sct status (CMD-05)
```

### Pattern 1: cli.py -> cli/__init__.py Migration (D-01, D-02)
**What:** Move existing `cli.py` contents to `cli/__init__.py`. Each submodule defines commands on its group's Typer app, imported and registered in `__init__.py`.
**When to use:** This is the foundational change -- must be done first.
**Example:**
```python
# cli/__init__.py (after migration)
from sc_tools.cli.qc import qc_app
from sc_tools.cli.preprocess import preprocess_app
from sc_tools.cli.validate import validate_app
from sc_tools.cli.benchmark import benchmark_app
from sc_tools.cli.status import status_app

app.add_typer(qc_app, name="qc")
app.add_typer(preprocess_app, name="preprocess")
# ...etc
```

```python
# cli/qc.py
import typer
from sc_tools.cli import cli_handler

qc_app = typer.Typer(help="Quality control commands")

@qc_app.callback(invoke_without_command=True)
def qc_callback(ctx: typer.Context) -> None:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())

@qc_app.command("run")
@cli_handler
def qc_run(
    file: str = typer.Argument(..., help="Path to h5ad checkpoint"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory"),
    modality: str = typer.Option("visium", help="Data modality"),
    sample_col: str = typer.Option("library_id", help="Sample column name"),
) -> CLIResult:
    # Lazy imports here
    ...
```

**Critical detail for D-02:** The old import path `from sc_tools.cli import app` must continue to work. Python resolves `sc_tools.cli` to `sc_tools/cli/__init__.py` automatically when converting a module to a package, so this is a natural consequence of the migration. However, any code that did `import sc_tools.cli` and accessed module-level attributes (e.g., `sc_tools.cli.qc_app`) will still work since `__init__.py` re-exports them.

### Pattern 2: Command Implementation Pattern
**What:** Every command follows the same structure: Typer args -> dep check -> lazy import -> load data -> call backend -> build CLIResult.
**Example:**
```python
@qc_app.command("run")
@cli_handler
def qc_run(
    file: str = typer.Argument(...),
    project_dir: str = typer.Option(".", "--project-dir"),
) -> CLIResult:
    from pathlib import Path
    _check_deps(["scanpy"])  # CMD-08: fast-fail

    import scanpy as sc
    from sc_tools.qc import compute_sample_metrics, classify_samples

    path = Path(file)
    if not path.exists():
        raise SCToolsUserError(f"File not found: {path}", suggestion="Check file path")

    adata = sc.read_h5ad(path)
    # ... call backend ...

    return CLIResult(
        status=Status.success,
        command="qc run",
        data={"metrics": metrics_dict, "classified": classified_dict},
        artifacts=[str(output_path)],
        provenance=Provenance(command="qc run"),
        message=f"QC metrics computed for {adata.n_obs} cells across {n_samples} samples",
    )
```

### Pattern 3: Fast-Fail Dependency Check (CMD-08)
**What:** A shared utility function that checks for optional dependencies before loading data.
**Example:**
```python
# In cli/__init__.py
_DEP_INSTALL = {
    "scanpy": "pip install scanpy",
    "scvi-tools": "pip install scvi-tools",
    "scib-metrics": "pip install scib-metrics",
    "h5py": "pip install h5py",
    "rapids_singlecell": "pip install rapids-singlecell (requires CUDA)",
}

def _check_deps(deps: list[str]) -> None:
    """Check optional dependencies before loading data. Raises SCToolsUserError."""
    missing = []
    for dep in deps:
        try:
            __import__(dep.replace("-", "_"))
        except ImportError:
            install = _DEP_INSTALL.get(dep, f"pip install {dep}")
            missing.append(f"  {dep}: {install}")
    if missing:
        raise SCToolsUserError(
            f"Missing required dependencies:\n" + "\n".join(missing),
            suggestion="Install missing packages and retry",
        )
```

### Pattern 4: Modality Auto-Detection (D-08)
**What:** Determine modality from `adata.uns` metadata or panel size heuristic.
**Example:**
```python
def _detect_modality(adata) -> str:
    """Auto-detect modality from adata.uns or panel size."""
    # Check adata.uns for explicit modality
    for key in ("modality", "platform", "data_type"):
        if key in adata.uns:
            val = str(adata.uns[key]).lower()
            if val in VALID_MODALITIES:
                return val

    # Heuristic: small panel -> targeted (xenium/cosmx)
    if adata.n_vars < 1000:
        logger.info("Auto-detected targeted panel (n_vars=%d < 1000)", adata.n_vars)
        return "xenium"  # generic targeted panel default

    return "visium"  # default for large panels
```

### Pattern 5: Shared Backend for CLI + MCP (CMD-07)
**What:** Both CLI commands and MCP tools call the same backend function, each wrapping the CLIResult differently.
**Example:**
```python
# Shared backend (e.g., in cli/validate.py or a separate module)
def _validate_backend(uri: str, phase: str, fix: bool = False) -> CLIResult:
    from sc_tools.validate import validate_file
    issues = validate_file(uri, phase=phase, fix=fix)
    return CLIResult(
        status=Status.success if not issues else Status.error,
        command=f"validate {phase}",
        data={"issues": issues, "phase": phase, "uri": uri},
        provenance=Provenance(command=f"validate {phase}"),
        message=f"Validation {'passed' if not issues else 'failed'}: {len(issues)} issues",
    )

# CLI: @cli_handler wraps and emits JSON to stdout
@validate_app.command()
@cli_handler
def validate_cmd(phase: str, file: str) -> CLIResult:
    return _validate_backend(file, phase)

# MCP: calls same backend, returns JSON string
@mcp.tool()
def validate_checkpoint(uri: str, phase: str, fix: bool = False) -> str:
    result = _validate_backend(uri, phase, fix)
    return result.model_dump_json()
```

### Pattern 6: sct status Without Registry (D-09)
**What:** `sct status` shows the phase DAG from `pipeline.py` even without a registry. With registry, shows completed phases and available next.
**Example:**
```python
@status_app.command()
@cli_handler
def status(project_dir: str = typer.Option(".", "--project-dir")) -> CLIResult:
    from sc_tools.pipeline import get_dag, get_available_next, get_phase_checkpoint

    dag = get_dag()
    completed = []
    registry_available = False

    try:
        # Try to read from registry
        # ... registry query ...
        registry_available = True
    except Exception:
        pass  # Graceful fallback

    available_next = get_available_next(completed)

    data = {
        "phases": {slug: {"label": spec.label, "checkpoint": spec.checkpoint}
                   for (_, slug), spec in dag.items()},
        "completed": [f"{g}/{s}" for g, s in completed],
        "available_next": [f"{g}/{s}" for g, s in available_next],
        "registry_available": registry_available,
    }

    return CLIResult(
        status=Status.success,
        command="status",
        data=data,
        provenance=Provenance(command="status"),
        message=f"{len(completed)} phases complete, {len(available_next)} available"
              + ("" if registry_available else " (registry unavailable)"),
    )
```

### Anti-Patterns to Avoid
- **Top-level heavy imports in cli/ submodules:** All scanpy, torch, scvi imports must be inside command functions (CLI-08)
- **Circular imports:** cli submodules import from `sc_tools.cli` (the `__init__`) for `cli_handler` and helpers. Do NOT have `__init__` import command implementations at module level -- use lazy registration or import inside the registration block
- **Direct print() in commands:** All output goes through CLIResult -> `_emit()`. Logging goes to stderr via Python `logging` module
- **Interactive prompts:** Violates CLI-07. All params via flags, fail fast on missing required args

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| QC metrics computation | Custom metric calculation | `compute_sample_metrics()` + `classify_samples()` | Handles modality defaults, MAD-based outliers, 15+ metrics |
| Preprocessing dispatch | Custom normalization pipeline | `preprocess()` from `recipes.py` | Handles 6 modalities, 4 integration methods, strategy pattern |
| Integration benchmarking | Custom embedding comparison | `compare_integrations()` from `integration.py` | Handles NaN masking, stratified subsampling, scib metrics |
| Checkpoint validation | Custom h5ad inspection | `validate_file()` from `validate.py` | Handles all 5 phases, legacy phase codes, auto-fix |
| Pipeline DAG traversal | Custom dependency tracking | `get_dag()` + `get_available_next()` from `pipeline.py` | Handles iterative phases, custom extensions, backward compat |
| HTML report generation | Custom HTML | `generate_{type}_report()` from `qc/report.py` and `bm/report.py` | Handles Bootstrap templates, Plotly embedding, tab comparison |
| Arrow string coercion | Custom dtype handling | Port `_coerce_arrow_strings()` from existing scripts | Prevents h5py serialization errors with newer pandas |

**Key insight:** Every CLI command is a thin wrapper. The backend functions are mature, tested, and handle edge cases. The CLI layer's job is: parse args, check deps, load data, call backend, package CLIResult.

## Common Pitfalls

### Pitfall 1: Circular Imports in cli/ Package
**What goes wrong:** `cli/__init__.py` imports from submodules which import `cli_handler` from `cli/__init__.py`.
**Why it happens:** Python package initialization order.
**How to avoid:** Have submodules import `cli_handler` and helpers from `sc_tools.cli` (absolute import). Register sub-apps in `__init__.py` after all definitions. The import chain is: submodule imports from `sc_tools.cli` -> `__init__.py` defines `cli_handler` before importing submodules.
**Warning signs:** `ImportError` on `from sc_tools.cli import app` or `from sc_tools.cli import cli_handler`.

### Pitfall 2: Arrow String Coercion
**What goes wrong:** AnnData loaded from h5ad on newer pandas gets ArrowStringArray columns, which h5py cannot serialize when writing back.
**Why it happens:** pandas >= 3.0 defaults to Arrow-backed strings.
**How to avoid:** Every command that loads h5ad should call `_coerce_arrow_strings()` (copy pattern from existing scripts). Set `pd.set_option("future.infer_string", False)` at module level in data-loading commands.
**Warning signs:** `TypeError` or `ArrowInvalid` when writing h5ad.

### Pitfall 3: CliRunner vs SystemExit in Tests
**What goes wrong:** `cli_handler` raises `SystemExit(0)` on success, which Typer's `CliRunner` catches and reports as exit_code=0. But `SystemExit(1/2/3)` for errors may cause test confusion.
**Why it happens:** `cli_handler` deliberately exits to ensure clean process termination.
**How to avoid:** In tests, check `result.exit_code` and parse `result.output` (stdout) as JSON to verify CLIResult contents. Use `CliRunner(mix_stderr=False)` to separate stderr from stdout.
**Warning signs:** Tests see non-zero exit codes when testing error paths.

### Pitfall 4: Typer Subcommand Naming with add_typer
**What goes wrong:** When using `app.add_typer(sub_app, name="qc")`, the subcommand `@sub_app.command("run")` becomes `sct qc run`. But if callback is misconfigured, `sct qc` alone may error instead of showing help.
**Why it happens:** `invoke_without_command=True` must be on the callback.
**How to avoid:** Keep the existing callback pattern from Phase 2 (`@qc_app.callback(invoke_without_command=True)`).
**Warning signs:** `sct qc` without subcommand exits with error instead of showing help.

### Pitfall 5: Benchmark Subsampling Transparency (D-12)
**What goes wrong:** Benchmark results don't clearly communicate whether/how subsampling was applied, leading to confusion when comparing across runs.
**Why it happens:** Subsampling decision happens deep in `compare_integrations()` but needs to surface in three places.
**How to avoid:** After calling `compare_integrations()`, inspect the returned DataFrame and provenance to record subsampling params. Include in CLIResult.provenance (extend Provenance model or use data dict) and pass to report generation.
**Warning signs:** User cannot determine from output whether results were subsampled.

### Pitfall 6: Large File Loading for QC
**What goes wrong:** `sct qc run` loads full AnnData into memory for datasets with millions of cells.
**Why it happens:** `sc.read_h5ad()` loads everything.
**How to avoid:** For Phase 3, accept the full load (Phase 7 MEM-01 addresses IO Gateway). Log cell count on load so user knows what's happening.
**Warning signs:** OOM on large datasets. `cli_handler` catches MemoryError and suggests `--subsample-n`.

## Code Examples

### cli/ Package __init__.py Structure
```python
# sc_tools/cli/__init__.py
"""Typer CLI application for sc_tools."""
from __future__ import annotations
import functools
import sys
import typer
from sc_tools.errors import (
    SCToolsDataError, SCToolsFatalError, SCToolsRuntimeError, SCToolsUserError,
)
from sc_tools.models.result import CLIResult, ErrorInfo, Provenance, Status

app = typer.Typer(name="sct", help="sct -- sc_tools command-line interface.",
                  pretty_exceptions_enable=False, add_completion=False)

_state: dict[str, bool] = {"human": False}

# ... _version_callback, main callback, _emit, _render_rich, cli_handler ...
# ... _check_deps, _make_error_result ...

# --- Register subcommands (after all helpers are defined) ---
from sc_tools.cli.qc import qc_app          # noqa: E402
from sc_tools.cli.preprocess import preprocess_app  # noqa: E402
from sc_tools.cli.validate import validate_app      # noqa: E402
from sc_tools.cli.benchmark import benchmark_app    # noqa: E402
from sc_tools.cli.status import status_app          # noqa: E402

app.add_typer(qc_app, name="qc")
app.add_typer(preprocess_app, name="preprocess")
# integrate_app stays as stub (no commands in Phase 3)
integrate_app = typer.Typer(help="Integration commands")
@integrate_app.callback(invoke_without_command=True)
def integrate_callback(ctx: typer.Context) -> None:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
app.add_typer(integrate_app, name="integrate")
app.add_typer(benchmark_app, name="benchmark")
# celltype stays as stub
celltype_app = typer.Typer(help="Cell typing commands")
@celltype_app.callback(invoke_without_command=True)
def celltype_callback(ctx: typer.Context) -> None:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
app.add_typer(celltype_app, name="celltype")
```

### sct benchmark integration Command
```python
# sc_tools/cli/benchmark.py
@benchmark_app.command("integration")
@cli_handler
def benchmark_integration(
    from_dir: str = typer.Option(..., "--from-dir", help="Directory with per-method h5ad files"),
    project_dir: str = typer.Option(".", "--project-dir"),
    batch_key: str = typer.Option("sample", "--batch-key"),
    celltype_key: str = typer.Option(None, "--celltype-key"),
    subsample_n: int = typer.Option(500_000, "--subsample-n"),
    seed: int = typer.Option(0, "--seed"),
    batch_weight: float = typer.Option(0.4, "--batch-weight"),
    bio_weight: float = typer.Option(0.6, "--bio-weight"),
) -> CLIResult:
    _check_deps(["h5py", "scanpy"])

    from pathlib import Path
    from sc_tools.bm.integration import compare_integrations
    # ... discover h5ad files in from_dir ...
    # ... build embedding_files dict ...
    # ... log subsampling if applied ...
    # ... call compare_integrations(embedding_files=..., subsample_n=..., seed=...) ...
    # ... build CLIResult with metrics, record subsampling in provenance ...
```

### Shared Test Fixture (conftest.py)
```python
# sc_tools/tests/conftest.py
import numpy as np
import pytest
from anndata import AnnData

@pytest.fixture
def adata_100():
    """100-cell AnnData for CLI integration tests (TST-05)."""
    np.random.seed(42)
    n_obs, n_vars = 100, 200
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype("float32")
    var_names = [f"MT-{i}" for i in range(3)] + [f"GENE{i}" for i in range(3, n_vars)]
    adata = AnnData(X)
    adata.var_names = var_names
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.obs["library_id"] = [f"L{i % 2}" for i in range(n_obs)]
    adata.obs["sample"] = adata.obs["library_id"]
    adata.obs["raw_data_dir"] = [f"/data/L{i % 2}" for i in range(n_obs)]
    adata.obsm["spatial"] = np.random.rand(n_obs, 2)
    return adata
```

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest >= 7.0 (installed) |
| Config file | `pyproject.toml` `[tool.pytest.ini_options]` |
| Quick run command | `pytest sc_tools/tests/test_cli.py -x -q` |
| Full suite command | `pytest sc_tools/tests/ -x -q` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| CMD-01 | `sct qc run` produces JSON metrics | integration | `pytest sc_tools/tests/test_cli_commands.py::test_qc_run -x` | Wave 0 |
| CMD-02 | `sct preprocess run` dispatches recipe | integration | `pytest sc_tools/tests/test_cli_commands.py::test_preprocess_run -x` | Wave 0 |
| CMD-03 | `sct validate` returns pass/fail | integration | `pytest sc_tools/tests/test_cli_commands.py::test_validate -x` | Wave 0 |
| CMD-04 | `sct benchmark integration` outputs ranked JSON | integration | `pytest sc_tools/tests/test_cli_commands.py::test_benchmark_integration -x` | Wave 0 |
| CMD-05 | `sct status` shows DAG state | integration | `pytest sc_tools/tests/test_cli_commands.py::test_status -x` | Wave 0 |
| CMD-06 | `sct report` generates HTML | integration | `pytest sc_tools/tests/test_cli_commands.py::test_report -x` | Wave 0 |
| CMD-07 | CLI and MCP share Result type | unit | `pytest sc_tools/tests/test_cli_commands.py::test_shared_result_type -x` | Wave 0 |
| CMD-08 | Missing dep fast-fail | unit | `pytest sc_tools/tests/test_cli_commands.py::test_dep_check -x` | Wave 0 |
| TST-05 | CLI integration tests with fixtures | integration | `pytest sc_tools/tests/test_cli_commands.py -x` | Wave 0 |
| TST-06 | E2E with real data (HPC) | e2e | `pytest sc_tools/tests/test_cli_e2e.py -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_cli.py sc_tools/tests/test_cli_commands.py -x -q`
- **Per wave merge:** `pytest sc_tools/tests/ -x -q`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/conftest.py` -- shared 100-cell AnnData fixture for CLI integration tests
- [ ] `sc_tools/tests/test_cli_commands.py` -- CLI command integration tests (CMD-01 through CMD-08)
- [ ] `sc_tools/tests/test_cli_e2e.py` -- E2E tests with real data, `skipif` guarded (TST-06)

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `scripts/run_*.py` with argparse | `sct` CLI with Typer + CLIResult | Phase 2-3 (2026-03) | Scripts become redundant, structured JSON output |
| MCP tools with custom JSON construction | MCP tools sharing CLIResult backend | Phase 2-3 (2026-03) | Single implementation, consistent output format |
| Manual `python scripts/run_qc_report.py --adata ...` | `sct qc run file.h5ad` / `sct report pre_filter --adata file.h5ad` | Phase 3 (2026-03) | Agent-native, no throwaway scripts |

## Open Questions

1. **How should `sct report` dispatch to report types?**
   - What we know: Four report types exist (`pre_filter`, `post_filter`, `post_integration`, `post_celltyping`). Each has different required inputs (e.g., `post_filter` needs both pre and post AnnData).
   - What's unclear: Whether `sct report pre_filter --adata file.h5ad` or `sct report --type pre_filter --adata file.h5ad`.
   - Recommendation: Use `sct report <type>` as a positional argument (matches `sct validate <phase>` pattern). Each type can have different optional flags.

2. **Where do CLI and MCP shared backends live?**
   - What we know: D-10 says they share backend functions. MCP already wraps `validate_file()` directly.
   - What's unclear: Whether shared backends should be in `cli/*.py` or in separate modules.
   - Recommendation: Backend functions already exist in `sc_tools.qc`, `sc_tools.pp`, `sc_tools.bm`, `sc_tools.validate`. CLI commands and MCP tools both call these. No new shared layer needed -- just ensure MCP tools use the same function calls as CLI commands.

3. **How to handle Arrow string coercion consistently?**
   - What we know: Three scripts have copy-pasted `_coerce_arrow_strings()`. It must run before any h5ad write.
   - Recommendation: Move to a shared utility (e.g., `sc_tools.ingest.loaders` or `sc_tools.utils`). Import in each CLI command that loads h5ad. This is a Claude's discretion item.

## Sources

### Primary (HIGH confidence)
- `sc_tools/cli.py` -- Phase 2 CLI foundation (directly read)
- `sc_tools/models/result.py` -- CLIResult model (directly read)
- `sc_tools/errors.py` -- Exception hierarchy (directly read)
- `sc_tools/validate.py` -- Checkpoint validation (directly read)
- `sc_tools/pipeline.py` -- Phase DAG (directly read)
- `sc_tools/pp/recipes.py` -- Preprocessing dispatch (directly read)
- `sc_tools/bm/integration.py` -- Integration benchmark (directly read)
- `sc_tools/bm/report.py` -- Benchmark report generation (directly read)
- `sc_tools/qc/report.py` -- QC report generation (directly read)
- `sc_tools/qc/sample_qc.py` -- QC metrics computation (directly read)
- `sc_tools/mcp/tools_server.py` -- MCP tool patterns (directly read)
- `scripts/run_qc_report.py` -- QC report script (directly read)
- `scripts/run_preprocessing.py` -- Preprocessing script (directly read)
- `scripts/run_integration_benchmark.py` -- Benchmark script (directly read)
- `sc_tools/tests/test_cli.py` -- Existing CLI tests (directly read)

### Secondary (MEDIUM confidence)
- Typer documentation for sub-app patterns (verified against codebase usage)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already in use, no new deps needed
- Architecture: HIGH -- cli/ package pattern is standard Python, existing code provides all patterns
- Pitfalls: HIGH -- identified from reading actual codebase, existing scripts show known gotchas (Arrow coercion, etc.)

**Research date:** 2026-03-21
**Valid until:** 2026-04-21 (stable -- no external dependency changes expected)
