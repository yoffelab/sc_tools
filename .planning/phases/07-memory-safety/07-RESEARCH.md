# Phase 7: Memory Safety - Research

**Researched:** 2026-03-24
**Domain:** Memory-safe IO for large h5ad files, pre-execution estimation, dry-run validation
**Confidence:** HIGH

## Summary

Phase 7 adds memory safety to the `sct` CLI through three mechanisms: (1) an IO Gateway (`sc_tools/io/`) with tiered loading that avoids full AnnData loads for metadata and summary operations, (2) `sct estimate` for pre-execution memory/runtime prediction from h5 metadata, and (3) `--dry-run` for all data-touching commands. The codebase already has proven patterns for all three -- h5py metadata reads in `inspect_checkpoint()`, h5py selective loads in `_load_embedding_h5py()`, backed mode in `smart_read_h5ad()`, and memory profiling in `memory/profiling.py`. The primary engineering work is consolidating these into the IO Gateway and wiring the CLI.

The key insight is that h5ad files are HDF5, so h5py can read shape, dtype, column names, and obsm keys without loading any array data. A 25GB file's metadata reads in milliseconds. Memory estimation before load is straightforward: `n_obs * n_vars * dtype.itemsize` for dense X, or inspecting sparse chunk sizes for CSR/CSC. The existing `estimate_adata_memory()` function operates on loaded AnnData -- the Gateway needs an equivalent that works on h5py file handles.

**Primary recommendation:** Build the IO Gateway as a thin wrapper around h5py with three entry points (`read_metadata`, `read_summary`, `read_full`) corresponding to tiers T1/T2/T3. Wire `cli_handler` to intercept `--dry-run` globally and check memory guard before T3 loads.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- D-01: New `sc_tools/io/` package -- dedicated module with gateway logic. Separate from `storage.py` (which handles fsspec URI resolution). IO Gateway focuses on memory-aware loading strategy.
- D-02: CLI-path only -- Gateway is used by `sct` CLI commands (via `cli_handler`). Library-level `sc_tools` functions keep calling `sc.read_h5ad` directly. Minimizes blast radius; agents get memory safety through CLI.
- D-03: Operation-based tier dispatch -- each CLI command declares what access level it needs: T1 (h5py metadata-only), T2 (backed/lazy mode), T3 (full load).
- D-04: Gateway reads the command's tier declaration and loads accordingly. No size-based auto-promotion.
- D-05: Pre-load memory guard with `--force` override. Gateway estimates memory from h5 metadata before T3 full load. If estimated memory exceeds available RAM, refuse with actionable error message suggesting `--force` flag or backed-mode alternative.

### Claude's Discretion
- How `sct estimate` calculates memory/runtime projections (formula-based from cell x gene counts, method-specific profiles, or empirical lookup)
- Internal structure of tier declaration (decorator parameter, enum, or command metadata)
- How `--dry-run` is implemented (global flag on `cli_handler` intercepting before execution, or per-command logic)
- What `--dry-run` reports (input validation, planned operations, estimated resources, or all of the above)
- How backed mode (T2) interacts with operations that need partial X access (e.g., QC metrics on a subset of genes)
- Whether to reuse/extend `inspect_checkpoint()` h5py pattern from MCP or write fresh gateway code
- Memory estimation formula calibration (static multipliers vs empirical profiling)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| MEM-01 | IO Gateway with tiered loading strategy: h5py for metadata/embeddings, backed mode for summaries, full load only for compute | Gateway architecture pattern, h5py metadata reading, backed mode via anndata, tier enum + decorator, existing codebase patterns from `inspect_checkpoint()` and `_load_embedding_h5py()` |
| MEM-02 | `sct estimate <command> <args>` -- pre-execution estimation of peak memory and runtime based on cell/gene count and method | h5py metadata extraction for n_obs/n_vars/dtype, memory formula (dense/sparse), method-specific multipliers, runtime heuristics |
| MEM-03 | `--dry-run` flag for all data-touching commands -- validate inputs, report what would happen, without executing | `cli_handler` global interception pattern, dry-run CLIResult envelope, existing dry_run patterns in `bm/slurm.py` |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| h5py | 3.15.1 | Low-level HDF5 access for T1 metadata reads | Already installed, proven in `_load_embedding_h5py` and `inspect_checkpoint` |
| anndata | (installed) | Backed-mode T2 reads via `read_h5ad(backed="r")` | Already used in `smart_read_h5ad`, native backed mode support |
| psutil | 7.2.2 | System memory detection for guard threshold | Already installed, used in `memory/profiling.py` |
| typer | 0.24.1 | CLI framework, global option injection for `--dry-run` | Already the CLI framework |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| scipy.sparse | (installed) | CSR/CSC matrix handling for sparse memory estimation | When X is stored as sparse in h5ad |
| numpy | (installed) | dtype introspection for memory calculation | Always, for `np.dtype(x).itemsize` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| h5py direct | zarr | zarr would add dep; h5ad files are HDF5 natively |
| psutil for RAM | `/proc/meminfo` | psutil is cross-platform (macOS + Linux), already installed |
| Backed mode (T2) | dask-backed AnnData | Complexity not justified; backed mode is sufficient for T2 use cases |

**Installation:**
```bash
# No new dependencies -- all are already installed
```

## Architecture Patterns

### Recommended Project Structure
```
sc_tools/
├── io/
│   ├── __init__.py          # Public API: read_metadata, read_summary, read_full
│   ├── gateway.py           # IOGateway class: tier dispatch, memory guard
│   ├── metadata.py          # T1: h5py metadata extraction (shape, cols, keys, dtype)
│   └── estimate.py          # Memory/runtime estimation from h5 metadata
├── cli/
│   ├── __init__.py          # cli_handler gains --dry-run + memory guard
│   ├── estimate.py          # NEW: sct estimate <command> <args>
│   └── [qc|preprocess|...].py  # Add tier declaration to each command
└── memory/
    └── profiling.py         # Existing -- unchanged
```

### Pattern 1: Tier Enum + Command Metadata

**What:** Each CLI command declares its data access tier via a parameter on `@cli_handler`.
**When to use:** Every data-touching CLI command.
**Example:**
```python
from enum import Enum

class DataTier(str, Enum):
    """Data access tier for IO Gateway (D-03)."""
    T1_METADATA = "metadata"   # h5py: shape, column names, obsm keys
    T2_SUMMARY = "summary"     # backed mode: obs/var DataFrames, obsm arrays
    T3_FULL = "full"           # full AnnData in memory

# In cli_handler decorator:
def cli_handler(func=None, *, tier: DataTier = DataTier.T3_FULL):
    """Wrap CLI command with error handling, tier dispatch, memory guard, dry-run."""
    ...
```

**Recommendation:** Use a decorator parameter on `cli_handler` rather than a separate decorator. This keeps the interception point single and consistent. The tier defaults to T3 (full load) for backwards compatibility.

### Pattern 2: IO Gateway with Pre-load Guard

**What:** Gateway checks available memory before T3 loads and refuses with actionable error if insufficient.
**When to use:** Every T3 (full load) operation.
**Example:**
```python
class IOGateway:
    """Memory-aware IO dispatcher (D-01, D-05)."""

    def read(self, path: str, tier: DataTier, *, force: bool = False):
        if tier == DataTier.T1_METADATA:
            return self._read_metadata(path)
        elif tier == DataTier.T2_SUMMARY:
            return self._read_summary(path)
        elif tier == DataTier.T3_FULL:
            if not force:
                self._check_memory_guard(path)
            return self._read_full(path)

    def _check_memory_guard(self, path: str) -> None:
        meta = self._read_metadata(path)
        estimated_mb = meta["estimated_memory_mb"]
        available_mb = psutil.virtual_memory().available / (1024 * 1024)
        if estimated_mb > available_mb * 0.8:  # 80% threshold
            raise SCToolsRuntimeError(
                f"Estimated memory ({estimated_mb:.0f}MB) exceeds 80% of available RAM ({available_mb:.0f}MB)",
                suggestion="Use --force to override, or reduce dataset size",
            )
```

### Pattern 3: Pre-load Memory Estimation from h5py

**What:** Estimate memory from h5ad file metadata without loading any array data.
**When to use:** Memory guard checks and `sct estimate`.
**Example:**
```python
def estimate_memory_from_h5(path: str) -> dict:
    """Estimate peak memory from h5ad metadata (no data load)."""
    import h5py

    with h5py.File(path, "r") as f:
        # X matrix estimation
        x_group = f["X"]
        if isinstance(x_group, h5py.Dataset):
            # Dense matrix
            x_bytes = x_group.shape[0] * x_group.shape[1] * x_group.dtype.itemsize
        elif "data" in x_group:
            # Sparse (CSR/CSC): data + indices + indptr
            x_bytes = (
                x_group["data"].shape[0] * x_group["data"].dtype.itemsize
                + x_group["indices"].shape[0] * x_group["indices"].dtype.itemsize
                + x_group["indptr"].shape[0] * x_group["indptr"].dtype.itemsize
            )

        # obs/var estimation (heuristic: ~200 bytes per cell for obs DataFrame)
        n_obs = ... # from shape
        obs_bytes = n_obs * 200  # conservative per-cell estimate

        # obsm estimation
        obsm_bytes = 0
        if "obsm" in f:
            for key in f["obsm"]:
                ds = f["obsm"][key]
                obsm_bytes += ds.shape[0] * ds.shape[1] * ds.dtype.itemsize

    total_mb = (x_bytes + obs_bytes + obsm_bytes) / (1024 * 1024)
    # Apply overhead multiplier (pandas conversion, anndata internals)
    return {"estimated_mb": total_mb * 1.5, "x_mb": x_bytes / (1024**2), ...}
```

### Pattern 4: Global --dry-run via cli_handler

**What:** `--dry-run` flag intercepted in `cli_handler` before command execution.
**When to use:** All data-touching commands.
**Example:**
```python
def cli_handler(func=None, *, tier: DataTier = DataTier.T3_FULL):
    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            dry_run = kwargs.pop("dry_run", False)
            if dry_run:
                # Validate inputs exist
                # Read T1 metadata
                # Report planned operations
                # Return CLIResult with status="dry_run" (or "skipped")
                result = _build_dry_run_result(fn, kwargs, tier)
                _emit(result)
                raise SystemExit(0)
            # ... normal execution
        return wrapper
    ...
```

**Recommendation for --dry-run injection:** Add `dry_run: bool = typer.Option(False, "--dry-run", help="...")` to each data-touching command's signature. The `cli_handler` detects it via `kwargs.get("dry_run")` and intercepts before calling the function body. This keeps Typer's help text accurate per-command.

### Pattern 5: Estimate Command

**What:** `sct estimate <command> <args>` parses the same arguments as `sct <command> <args>` but only reads T1 metadata and returns projections.
**When to use:** Before running expensive commands on large datasets.
**Example:**
```python
# sct estimate preprocess run /path/to/big.h5ad --modality visium
# Returns JSON with estimated_peak_memory_mb, estimated_runtime_s
```

**Recommendation:** Implement as a top-level command that takes the target command name and its file argument, reads T1 metadata, and applies method-specific multipliers. It does NOT re-parse all args of the target command -- just needs the file path and optionally the method/command name for method-specific estimates.

### Anti-Patterns to Avoid
- **Loading full AnnData for metadata queries:** Use h5py directly. `sc.read_h5ad` for a 25G file allocates 25G+ of RAM.
- **Auto-promoting tiers based on file size:** D-04 explicitly rejects this. Tier is operation-based, not data-based.
- **Memory guard without override:** D-05 requires `--force` flag for expert users who know their system.
- **Modifying library-level sc_tools code:** D-02 restricts Gateway to CLI path only.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| HDF5 metadata reading | Custom binary parser | h5py | HDF5 format is complex; h5py is the standard Python binding |
| System memory detection | `/proc/meminfo` parser | psutil | Cross-platform (macOS + Linux), already installed |
| Backed-mode AnnData | Custom lazy loader | `anndata.read_h5ad(backed="r")` | Native support, handles sparse/dense/categorical automatically |
| Sparse matrix size estimation | Manual byte counting | Read h5py dataset shapes + dtypes | h5py exposes shape and dtype without loading data |

**Key insight:** The h5ad format is well-structured HDF5. All metadata needed for memory estimation is accessible via h5py attributes and dataset shapes without reading any array data. The existing codebase already has two proven h5py read patterns to consolidate.

## Common Pitfalls

### Pitfall 1: Categorical Columns in h5py
**What goes wrong:** obs columns stored as categorical have a `codes` + `categories` structure in HDF5, not plain arrays. Naive reads get raw integer codes.
**Why it happens:** AnnData serializes `pd.Categorical` as two datasets.
**How to avoid:** The `_load_embedding_h5py()` pattern at line 101-107 already handles this -- reuse that pattern in the Gateway.
**Warning signs:** Getting integer arrays where you expect string labels.

### Pitfall 2: Sparse X Format Detection
**What goes wrong:** Assuming X is always a dense dataset. Many h5ad files store X as CSR/CSC sparse.
**Why it happens:** anndata writes sparse matrices as HDF5 groups with `data`, `indices`, `indptr` sub-datasets rather than a single dataset.
**How to avoid:** Check if `f["X"]` is an `h5py.Group` (sparse) vs `h5py.Dataset` (dense). For sparse, sum component sizes.
**Warning signs:** `KeyError` when accessing `f["X"].shape` on a group, or memory estimates wildly wrong.

### Pitfall 3: Backed Mode Limitations
**What goes wrong:** Operations that need X values (e.g., `sc.pp.normalize_total`) fail on backed AnnData because they try to write to the read-only array.
**Why it happens:** Backed mode keeps data on disk. Write operations require explicit `.to_memory()` first.
**How to avoid:** T2 should only be used for reading obs/var/obsm. Commands needing X computation must use T3. This is exactly why tier is operation-based (D-03/D-04).
**Warning signs:** `ValueError: assignment destination is read-only` during backed-mode operations.

### Pitfall 4: Memory Estimation Accuracy
**What goes wrong:** Underestimating memory because the formula only counts X + obs + obsm, missing anndata overhead, pandas index duplication, and temporary allocations during scanpy operations.
**Why it happens:** Peak memory during preprocessing can be 2-3x the final AnnData size due to intermediate copies.
**How to avoid:** Apply a 1.5-2x overhead multiplier for T3 estimates. For preprocessing specifically (normalize + HVG + PCA + integration), use 2.5-3x because of intermediate dense copies.
**Warning signs:** OOM despite estimate saying it should fit.

### Pitfall 5: n_obs Detection in h5py
**What goes wrong:** Different h5ad file versions store shape information differently. Some have `_index_length` attribute, others use `_index` dataset length.
**Why it happens:** anndata format has evolved. The `inspect_checkpoint()` function at line 510-523 already handles three fallback strategies.
**How to avoid:** Reuse the fallback pattern from `inspect_checkpoint()`: attrs `_index` -> dataset `_index` length -> first key length.
**Warning signs:** Getting `None` or `?` for n_obs.

### Pitfall 6: --dry-run and Typer Callback Ordering
**What goes wrong:** Adding `--dry-run` as a global option on the app callback means it's not available in command-specific help.
**Why it happens:** Typer's callback options are separate from command options.
**How to avoid:** Add `--dry-run` as a parameter to each data-touching command. `cli_handler` intercepts it from kwargs.
**Warning signs:** `--dry-run` not showing up in `sct qc run --help`.

## Code Examples

### T1 Metadata Read (from inspect_checkpoint pattern)
```python
# Source: sc_tools/mcp/tools_server.py lines 509-536
def read_h5ad_metadata(path: str) -> dict:
    """Read metadata from h5ad without loading arrays."""
    import h5py

    with h5py.File(path, "r") as f:
        # Shape detection with fallbacks
        n_obs = f["obs"].attrs.get("_index_length", None)
        if n_obs is None and "_index" in f["obs"]:
            n_obs = len(f["obs"]["_index"])
        elif n_obs is None:
            first_key = next(iter(f["obs"].keys()), None)
            n_obs = len(f["obs"][first_key]) if first_key else 0

        n_vars = f["var"].attrs.get("_index_length", None)
        if n_vars is None and "_index" in f["var"]:
            n_vars = len(f["var"]["_index"])
        elif n_vars is None:
            first_key = next(iter(f["var"].keys()), None)
            n_vars = len(f["var"][first_key]) if first_key else 0

        # X dtype and sparsity
        x_node = f["X"]
        is_sparse = isinstance(x_node, h5py.Group)
        if is_sparse:
            x_dtype = str(x_node["data"].dtype)
        else:
            x_dtype = str(x_node.dtype)

        obs_cols = [k for k in f["obs"].keys() if not k.startswith("__")]
        obsm_keys = list(f["obsm"].keys()) if "obsm" in f else []
        layer_keys = list(f["layers"].keys()) if "layers" in f else []
        uns_keys = list(f["uns"].keys()) if "uns" in f else []

    return {
        "n_obs": n_obs, "n_vars": n_vars,
        "x_dtype": x_dtype, "x_sparse": is_sparse,
        "obs_columns": obs_cols, "obsm_keys": obsm_keys,
        "layer_keys": layer_keys, "uns_keys": uns_keys,
    }
```

### Memory Estimation Formula
```python
# Source: adapted from sc_tools/memory/profiling.py estimate_adata_memory
def estimate_from_metadata(meta: dict) -> dict:
    """Estimate memory from T1 metadata dict."""
    import numpy as np

    n_obs = meta["n_obs"]
    n_vars = meta["n_vars"]
    dtype_size = np.dtype(meta["x_dtype"]).itemsize

    if meta["x_sparse"]:
        # Sparse: estimate ~10% density for scRNA-seq
        density = 0.1
        nnz = int(n_obs * n_vars * density)
        x_bytes = nnz * dtype_size + nnz * 4 + (n_obs + 1) * 4  # data + indices + indptr
    else:
        x_bytes = n_obs * n_vars * dtype_size

    obs_bytes = n_obs * 200  # ~200 bytes/cell for obs DataFrame
    var_bytes = n_vars * 100  # ~100 bytes/gene for var DataFrame

    base_mb = (x_bytes + obs_bytes + var_bytes) / (1024 ** 2)

    return {
        "x_mb": x_bytes / (1024 ** 2),
        "obs_mb": obs_bytes / (1024 ** 2),
        "base_mb": base_mb,
        "estimated_peak_mb": base_mb * 2.0,  # 2x overhead for processing
    }
```

### Command Tier Declaration
```python
# Source: new pattern for Phase 7
@qc_app.command("run")
@cli_handler(tier=DataTier.T3_FULL)
def qc_run(
    file: str = typer.Argument(...),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs and report plan without executing"),
    force: bool = typer.Option(False, "--force", help="Override memory guard"),
    ...
):
    ...
```

### Tier-to-Command Mapping
```python
# Tier assignments per CONTEXT.md D-03:
TIER_MAP = {
    "status show": DataTier.T1_METADATA,       # Shape, column names
    "validate run": DataTier.T1_METADATA,       # Structural checks only
    "provenance show": DataTier.T1_METADATA,    # Read uns/sidecar
    "estimate": DataTier.T1_METADATA,           # Pre-execution sizing

    "qc run": DataTier.T2_SUMMARY,             # Needs obs columns (but not full X for basic QC)
    "benchmark integration": DataTier.T2_SUMMARY,  # Reads obsm embeddings via h5py (already does this)

    "preprocess run": DataTier.T3_FULL,         # Needs full X for normalization
    "report generate": DataTier.T3_FULL,        # Needs full AnnData for plots
    "de run": DataTier.T3_FULL,                 # Needs raw counts
}
```

**Note on qc_run tier:** QC currently does `sc.read_h5ad` then `calculate_qc_metrics`. The QC metrics computation (`scanpy.pp.calculate_qc_metrics`) needs X values. However, if QC metrics are already computed (stored in obs columns), qc_run could work at T2. For initial implementation, keep at T3 to avoid breaking existing behavior; optimize later if needed.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `sc.read_h5ad()` for everything | h5py for metadata, backed mode for summaries | Established pattern in codebase | Prevents OOM on 25G+ files |
| No pre-load check | Memory guard with `--force` override | New (this phase) | Catches OOM before 44GB allocation |
| Manual `--dry-run` per module | Global `--dry-run` on all CLI commands | New (this phase) | Consistent agent experience |

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest (installed) |
| Config file | pyproject.toml `[tool.pytest.ini_options]` |
| Quick run command | `python -m pytest sc_tools/tests/test_io_gateway.py -x` |
| Full suite command | `python -m pytest sc_tools/tests/ -v --tb=short` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| MEM-01a | T1 reads metadata via h5py without full load | unit | `pytest sc_tools/tests/test_io_gateway.py::test_t1_metadata_read -x` | Wave 0 |
| MEM-01b | T2 reads obs/var via backed mode | unit | `pytest sc_tools/tests/test_io_gateway.py::test_t2_backed_read -x` | Wave 0 |
| MEM-01c | T3 full load with memory guard | unit | `pytest sc_tools/tests/test_io_gateway.py::test_t3_memory_guard -x` | Wave 0 |
| MEM-01d | Memory guard blocks when estimated > available | unit | `pytest sc_tools/tests/test_io_gateway.py::test_memory_guard_blocks -x` | Wave 0 |
| MEM-01e | --force overrides memory guard | unit | `pytest sc_tools/tests/test_io_gateway.py::test_force_override -x` | Wave 0 |
| MEM-02a | Estimate returns peak memory from h5 metadata | unit | `pytest sc_tools/tests/test_estimate.py::test_estimate_memory -x` | Wave 0 |
| MEM-02b | Estimate handles sparse X | unit | `pytest sc_tools/tests/test_estimate.py::test_estimate_sparse -x` | Wave 0 |
| MEM-02c | sct estimate CLI command returns JSON | integration | `pytest sc_tools/tests/test_cli_estimate.py::test_estimate_cli -x` | Wave 0 |
| MEM-03a | --dry-run validates inputs and exits 0 | integration | `pytest sc_tools/tests/test_cli_dryrun.py::test_dryrun_validates_inputs -x` | Wave 0 |
| MEM-03b | --dry-run reports planned operations | integration | `pytest sc_tools/tests/test_cli_dryrun.py::test_dryrun_reports_plan -x` | Wave 0 |
| MEM-03c | --dry-run does not modify data | integration | `pytest sc_tools/tests/test_cli_dryrun.py::test_dryrun_no_side_effects -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `python -m pytest sc_tools/tests/test_io_gateway.py sc_tools/tests/test_estimate.py sc_tools/tests/test_cli_dryrun.py -x`
- **Per wave merge:** `python -m pytest sc_tools/tests/ -v --tb=short`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/test_io_gateway.py` -- covers MEM-01 (T1/T2/T3 reads, memory guard, force override)
- [ ] `sc_tools/tests/test_estimate.py` -- covers MEM-02 (memory estimation from h5 metadata, sparse/dense)
- [ ] `sc_tools/tests/test_cli_estimate.py` -- covers MEM-02 CLI integration
- [ ] `sc_tools/tests/test_cli_dryrun.py` -- covers MEM-03 (--dry-run on data-touching commands)

## Open Questions

1. **QC tier assignment**
   - What we know: `qc_run` currently loads full AnnData for `calculate_qc_metrics` which needs X. But if metrics are already computed, obs columns suffice.
   - What's unclear: Should qc_run dynamically choose tier based on whether metrics exist?
   - Recommendation: Start at T3 for qc_run. Add a T2 fast path later if needed. Keep scope manageable.

2. **Runtime estimation accuracy**
   - What we know: Memory estimation from file metadata is reliable (h5py shape/dtype). Runtime estimation is much harder -- depends on CPU/GPU, disk speed, integration method convergence.
   - What's unclear: How accurate can runtime estimates be?
   - Recommendation: Return memory estimates with HIGH confidence and runtime estimates as "rough order of magnitude" with explicit LOW confidence label in the JSON output. Use simple formulas: `n_obs * n_vars * method_factor / reference_throughput`.

3. **Sparse density estimation**
   - What we know: h5py can read the actual nnz from `f["X/data"].shape[0]` for sparse matrices without loading data.
   - What's unclear: Whether all h5ad files store sparse X with consistent group structure.
   - Recommendation: Read actual nnz from h5py when available, fall back to 10% density estimate. The h5py approach is direct and accurate.

## Project Constraints (from CLAUDE.md)

- Output paths: active projects at `~/Documents/projects/active/<project>/` -- never repo root
- No heavy imports at module level (CLI-08) -- IO Gateway must lazy-import h5py
- `cli_handler` is the single interception point for error handling, provenance, and now dry-run/memory guard
- CLIResult envelope is the standard output format -- dry-run results use it too
- Follow `_check_deps()` pattern for h5py dependency checking
- ruff UP017: Use `timezone.utc` with `# noqa: UP017` for Python 3.10 compat

## Sources

### Primary (HIGH confidence)
- `sc_tools/mcp/tools_server.py` lines 487-549 -- proven h5py metadata read pattern
- `sc_tools/bm/integration.py` lines 63-112 -- proven h5py selective embedding load
- `sc_tools/storage.py` lines 154-177 -- backed mode AnnData loading
- `sc_tools/memory/profiling.py` -- existing memory profiling utilities
- `sc_tools/cli/__init__.py` lines 223-276 -- cli_handler decorator pattern
- `sc_tools/models/result.py` -- CLIResult envelope model
- `sc_tools/bm/slurm.py` lines 149-178 -- existing dry_run pattern

### Secondary (MEDIUM confidence)
- h5py 3.15.1 API: Dataset.shape, Dataset.dtype, Group vs Dataset detection
- psutil 7.2.2: `virtual_memory().available` for RAM detection

### Tertiary (LOW confidence)
- Runtime estimation multipliers (method-specific throughput factors) -- will need empirical calibration

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all dependencies already installed and proven in codebase
- Architecture: HIGH - consolidating existing proven patterns, decisions locked in CONTEXT.md
- Pitfalls: HIGH - identified from actual codebase patterns and h5ad format experience

**Research date:** 2026-03-24
**Valid until:** 2026-04-24 (stable domain, no external API changes expected)
