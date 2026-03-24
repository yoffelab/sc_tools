# Phase 5: Provenance - Research

**Researched:** 2026-03-23
**Domain:** CLI provenance tracking, file checksumming, lineage tracing
**Confidence:** HIGH

## Summary

Phase 5 expands the minimal `Provenance` Pydantic model (command, timestamp, version) into a full provenance system that auto-writes `.provenance.json` sidecars alongside every CLI artifact, embeds provenance in `adata.uns` for h5ad outputs, provides two new CLI commands (`sct provenance show` and `sct provenance trace`), and threads `random_state` through all Leiden clustering call paths for reproducibility.

The implementation is entirely within the existing Python stdlib and Pydantic stack -- no new external dependencies are needed. SHA256 checksumming uses `hashlib`, peak memory uses `resource.getrusage` (with platform-aware unit handling for macOS vs Linux), and timing uses `time.monotonic`. The `cli_handler` decorator is the single injection point for automatic sidecar writing after successful command execution.

**Primary recommendation:** Expand the `Provenance` model with nested Pydantic sub-models (`InputFile`, full `Provenance`), hook sidecar writing into `cli_handler` after `_emit()`, and register a new `provenance_app` Typer group following the `register_discovery(app)` pattern from Phase 4.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Auto sidecar writing via `cli_handler`. After successful execution, if `result.artifacts` is non-empty, write a `.provenance.json` sidecar next to each artifact automatically. No opt-out flag.
- **D-02:** Sidecar naming: `{original_filename}.provenance.json`.
- **D-03:** Only written on success -- error results do not generate sidecars.
- **D-04:** Hybrid approach: always write sidecar files; additionally embed in `adata.uns['sct_provenance']` for h5ad outputs.
- **D-05:** Sidecars are canonical for `sct provenance trace`. `adata.uns` is portable backup -- trace falls back to it when sidecar is missing.
- **D-06:** Expand existing `Provenance` model with: command, params, inputs (with SHA256), sc_tools_version, timestamp, runtime_s, peak_memory_mb.
- **D-07:** All CLI flags serialized as-passed including defaults.
- **D-08:** Input file records: path (relative to project root), path_type ("relative"), sha256, size_bytes.
- **D-09:** Store paths relative to project root (using `--project-dir` or `.` default).
- **D-10:** Missing input file resolution: check adata.uns, then SHA256 search in project dir.
- **D-11:** `sct provenance trace` follows sidecar input references recursively. Stop at raw data origin (no sidecar/uns).
- **D-12:** Trace output is flat chronological list with file path, command, timestamp, note field.
- **D-13:** `sct provenance show` reads and displays sidecar (or uns fallback) -- no recursion.
- **D-14:** Thread `random_state` through all `sc.tl.leiden()` call paths. Default `random_state=0`.
- **D-15:** Record resolution and random_state in provenance params. Test identical-params-identical-results.

### Claude's Discretion
- SHA256 computation strategy for large h5ad files (full file vs chunked streaming)
- Internal Provenance model structure (nested Pydantic models vs flat dict)
- How cli_handler passes input file info to provenance (CLIResult augmentation)
- How adata.uns['sct_provenance'] is structured
- Depth limit for trace recursion (or unlimited with cycle detection)
- How peak_memory_mb is measured

### Deferred Ideas (OUT OF SCOPE)
- W3C PROV-JSON export (INF-01)
- RO-Crate packaging (INF-02)
- Database-backed provenance (INF-03)
- DOT/graph visualization of lineage DAG
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| PRV-01 | JSON sidecar `.provenance.json` written alongside every CLI output file | cli_handler hook after _emit(), sidecar naming convention D-02, success-only guard D-03 |
| PRV-02 | Sidecar includes: command, params, input files with SHA256 checksums, sc_tools version, timestamp, runtime_s, peak_memory_mb | Expanded Provenance model with InputFile sub-model; resource.getrusage for memory; time.monotonic for timing |
| PRV-03 | `sct provenance show <file>` -- display provenance for single output | New provenance_app Typer group, register_provenance(app) pattern, sidecar read with adata.uns fallback |
| PRV-04 | `sct provenance trace <file>` -- trace full lineage DAG via input file references | Recursive sidecar walking, SHA256-based relocation, cycle detection via visited set |
| PRV-05 | Reproducible Leiden clustering -- configurable resolution and random_state | Thread random_state through _leiden_cluster, compute_integration_metrics, pp.reduce.leiden, pp.strategy |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pydantic | (existing) | Provenance model, InputFile model, serialization | Already used for CLIResult -- extend same pattern |
| hashlib | stdlib | SHA256 file checksumming | No external deps, streaming-capable |
| resource | stdlib | Peak memory measurement via getrusage | Cross-platform (macOS/Linux), no overhead |
| time | stdlib | Runtime measurement via monotonic clock | Wall-clock precision, no drift |
| pathlib | stdlib | Relative path computation | Already used throughout CLI |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| h5py | (existing) | Read adata.uns for provenance fallback in trace | Only for provenance show/trace on h5ad files |
| typer | (existing) | New provenance command group | Already the CLI framework |
| json | stdlib | Read/write sidecar files | Standard JSON I/O |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| resource.getrusage | tracemalloc | tracemalloc tracks Python allocations only, not C extensions (numpy, scanpy). getrusage captures total process RSS which is more accurate for bioinformatics workloads |
| Full-file SHA256 | Header-only hash | Header-only is faster for multi-GB h5ad but loses integrity guarantee. Full-file streaming (64KB chunks) is ~2 seconds for 1GB -- acceptable for CLI command that just processed the file |
| Nested Pydantic models | Flat dict | Nested models provide validation, serialization, and schema generation. Matches existing CLIResult pattern |

**Installation:** No new packages needed -- all stdlib or existing dependencies.

## Architecture Patterns

### Recommended Project Structure
```
sc_tools/
├── models/
│   └── result.py        # Expand Provenance, add InputFile, add ProvenanceRecord
├── cli/
│   ├── __init__.py      # Hook sidecar writing into cli_handler
│   └── provenance.py    # New: show and trace commands (register_provenance pattern)
├── provenance/          # New: provenance utilities module
│   ├── __init__.py
│   ├── sidecar.py       # write_sidecar(), read_sidecar(), sidecar_path_for()
│   ├── checksum.py      # sha256_file() streaming hasher
│   └── trace.py         # trace_lineage() recursive walker
├── bm/
│   └── integration.py   # Thread random_state through _leiden_cluster
└── pp/
    └── reduce.py        # Thread random_state through leiden()
```

### Pattern 1: Expanded Provenance Model (Nested Pydantic)
**What:** Replace the minimal Provenance model with a full ProvenanceRecord containing nested InputFile objects.
**When to use:** Every sidecar write, every provenance display.
**Example:**
```python
# In sc_tools/models/result.py
class InputFile(BaseModel):
    """Input file record for provenance tracking (D-08)."""
    path: str           # Relative to project root
    path_type: str = "relative"
    sha256: str
    size_bytes: int

class ProvenanceRecord(BaseModel):
    """Full provenance sidecar content (D-06)."""
    command: str
    params: dict[str, Any]
    inputs: list[InputFile] = Field(default_factory=list)
    sc_tools_version: str = Field(default_factory=_get_version)
    timestamp: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),  # noqa: UP017
    )
    runtime_s: float | None = None
    peak_memory_mb: float | None = None

# Keep existing Provenance class for backwards compat in CLIResult
# ProvenanceRecord is for sidecars only
```

### Pattern 2: cli_handler Sidecar Hook
**What:** After successful `_emit()`, if the CLIResult has artifacts, compute provenance and write sidecars.
**When to use:** Automatically for every @cli_handler-decorated command.
**Example:**
```python
def cli_handler(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.monotonic()
        mem_before = _get_peak_memory_mb()
        try:
            result = func(*args, **kwargs)
            runtime_s = time.monotonic() - start_time
            peak_mb = _get_peak_memory_mb()
            _emit(result)
            # Sidecar writing: only on success with artifacts (D-01, D-03)
            if result.status == Status.success and result.artifacts:
                _write_provenance_sidecars(result, kwargs, runtime_s, peak_mb)
            raise SystemExit(0)
        except SystemExit:
            raise
        # ... existing error handlers ...
    return wrapper
```

### Pattern 3: register_provenance(app) Registration
**What:** Follow the Phase 4 `register_discovery(app)` pattern to avoid circular imports.
**When to use:** Registering the provenance command group.
**Example:**
```python
# In sc_tools/cli/provenance.py
def register_provenance(target_app: typer.Typer) -> None:
    """Register provenance commands on the given Typer app."""
    provenance_app = typer.Typer(help="Provenance commands")

    @provenance_app.command("show")
    @cli_handler
    def provenance_show(file: str = typer.Argument(...)):
        ...

    @provenance_app.command("trace")
    @cli_handler
    def provenance_trace(file: str = typer.Argument(...)):
        ...

    target_app.add_typer(provenance_app, name="provenance")

# In sc_tools/cli/__init__.py (at bottom with other registrations)
from sc_tools.cli.provenance import register_provenance  # noqa: E402
register_provenance(app)
```

### Pattern 4: Input File Tracking via CLIResult Augmentation
**What:** Commands declare their input files in CLIResult.data so cli_handler can compute checksums.
**When to use:** Every data-producing command.
**Example:**
```python
# In a command function, include input files in data dict
return CLIResult(
    ...
    data={
        "input": str(path),         # existing pattern -- already used
        "_input_files": [str(path)],  # explicit list for provenance
        ...
    },
    artifacts=[str(output_path)],
    ...
)
```
**Recommended approach:** Use the existing `data["input"]` pattern plus a new convention key `_input_files` (list of paths). The cli_handler reads this key, computes SHA256 for each, builds InputFile records, and removes the key before emitting JSON. Alternatively, add an `input_files` field to CLIResult directly -- this is cleaner but requires updating all existing commands. The `_input_files` convention key is pragmatic for Phase 5 without touching Phase 2-4 code.

### Anti-Patterns to Avoid
- **Blocking on SHA256 for multi-GB files in the main thread:** Use streaming 64KB chunks; never read entire file into memory.
- **Absolute paths in sidecars:** Always compute relative paths using `--project-dir`. Absolute paths break when projects move between machines/scratch.
- **Writing sidecar before output file:** The artifact must exist before computing its companion sidecar. Write order: output file, then sidecar.
- **Modifying Provenance class incompatibly:** The existing `Provenance` in CLIResult must stay backwards compatible. Create a new `ProvenanceRecord` for sidecars.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SHA256 checksumming | Custom byte reading | `hashlib.sha256()` with streaming `file.read(65536)` | Edge cases: partial reads, encoding, huge files |
| Peak memory measurement | Custom `/proc/self/status` parser | `resource.getrusage(resource.RUSAGE_SELF).ru_maxrss` | Cross-platform, handles macOS (bytes) vs Linux (KB) automatically with platform check |
| JSON sidecar I/O | Custom file format | `json.dump/load` with Pydantic `.model_dump(mode="json")` | Schema validation, round-trip fidelity |
| Relative path computation | String manipulation | `pathlib.Path.relative_to()` | Handles edge cases (symlinks, `..` components) |

**Key insight:** The provenance system is entirely composed of well-understood stdlib operations. The complexity is in the integration points (cli_handler hook, input file tracking convention) not the individual operations.

## Common Pitfalls

### Pitfall 1: macOS vs Linux ru_maxrss Units
**What goes wrong:** `resource.getrusage().ru_maxrss` returns bytes on macOS but kilobytes on Linux. Dividing by 1024 on both platforms gives wrong results on one.
**Why it happens:** Different OS kernels, different conventions.
**How to avoid:**
```python
import platform, resource
ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if platform.system() == "Darwin":
    peak_mb = ru / (1024 * 1024)  # bytes -> MB
else:
    peak_mb = ru / 1024  # KB -> MB
```
**Warning signs:** peak_memory_mb values that are 1000x too large or small.

### Pitfall 2: Sidecar Race Condition with Parallel Execution
**What goes wrong:** If two CLI commands write to the same output path concurrently, sidecars can be corrupted.
**Why it happens:** JSON write is not atomic.
**How to avoid:** Write to a temp file then `os.replace()` (atomic on POSIX). This is a defense-in-depth measure -- the CLI is typically sequential.
```python
import tempfile, os
tmp_fd, tmp_path = tempfile.mkstemp(dir=sidecar_path.parent, suffix=".tmp")
with os.fdopen(tmp_fd, "w") as f:
    json.dump(record, f, indent=2)
os.replace(tmp_path, sidecar_path)
```
**Warning signs:** Truncated or mixed JSON in sidecar files.

### Pitfall 3: Circular Reference in Lineage Trace
**What goes wrong:** If a file's provenance lists itself as an input (or cycle in DAG), trace enters infinite recursion.
**Why it happens:** Unlikely in practice but possible with manual sidecar editing or bugs.
**How to avoid:** Maintain a `visited: set[str]` of canonical file paths (resolved via SHA256 or absolute path). Skip already-visited nodes.
**Warning signs:** `RecursionError` or extremely long trace output.

### Pitfall 4: Relative Path Resolution Failures
**What goes wrong:** `Path.relative_to()` raises ValueError when paths are not relative.
**Why it happens:** Artifact path not under project-dir (e.g., writing to /tmp or absolute path).
**How to avoid:** Use `try/except ValueError` and fall back to the absolute path. Log a warning.
**Warning signs:** `ValueError` in sidecar writing path.

### Pitfall 5: h5ad Read for adata.uns Fallback in Trace
**What goes wrong:** Loading full AnnData for a 10GB file just to read `adata.uns` provenance.
**Why it happens:** Using `sc.read_h5ad()` instead of targeted h5py access.
**How to avoid:** Use h5py to read only `uns/sct_provenance` key, not full AnnData.
```python
import h5py, json
with h5py.File(path, "r") as f:
    if "uns/sct_provenance" in f:
        raw = f["uns/sct_provenance"][()]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        prov = json.loads(raw)
```
**Warning signs:** Long loading times in `sct provenance trace`.

### Pitfall 6: scanpy.tl.leiden random_state Parameter
**What goes wrong:** `sc.tl.leiden()` accepts `random_state` but it may not guarantee determinism across all igraph/leidenalg versions.
**Why it happens:** Leiden algorithm depends on igraph's random number generator seeding.
**How to avoid:** Thread `random_state` and document that reproducibility is best-effort across different leidenalg versions. Test on the same environment (same leidenalg version).
**Warning signs:** Different cluster labels despite identical parameters on different machines.

## Code Examples

### SHA256 Streaming Checksum
```python
# In sc_tools/provenance/checksum.py
import hashlib
from pathlib import Path

def sha256_file(path: Path | str, chunk_size: int = 65536) -> str:
    """Compute SHA256 hex digest of a file using streaming reads."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()
```

### Peak Memory Helper
```python
# In sc_tools/provenance/sidecar.py
import platform
import resource

def get_peak_memory_mb() -> float:
    """Return peak RSS memory in megabytes."""
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() == "Darwin":
        return ru / (1024 * 1024)
    return ru / 1024  # Linux: KB -> MB
```

### Sidecar Writing
```python
# In sc_tools/provenance/sidecar.py
import json
import os
import tempfile
from pathlib import Path
from sc_tools.models.result import ProvenanceRecord

def sidecar_path_for(artifact_path: str | Path) -> Path:
    """Return the sidecar path for an artifact (D-02)."""
    return Path(str(artifact_path) + ".provenance.json")

def write_sidecar(artifact_path: str | Path, record: ProvenanceRecord) -> Path:
    """Write provenance sidecar atomically (D-01)."""
    sp = sidecar_path_for(artifact_path)
    data = record.model_dump(mode="json")
    # Atomic write
    tmp_fd, tmp_path = tempfile.mkstemp(dir=sp.parent, suffix=".tmp")
    try:
        with os.fdopen(tmp_fd, "w") as f:
            json.dump(data, f, indent=2)
        os.replace(tmp_path, sp)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise
    return sp

def read_sidecar(artifact_path: str | Path) -> dict | None:
    """Read provenance sidecar if it exists."""
    sp = sidecar_path_for(artifact_path)
    if not sp.exists():
        return None
    with open(sp) as f:
        return json.load(f)
```

### Embedding Provenance in adata.uns
```python
# After writing h5ad output, embed provenance (D-04)
import json

def embed_provenance_in_adata(adata, record: ProvenanceRecord) -> None:
    """Embed provenance in adata.uns for portability (D-04)."""
    adata.uns["sct_provenance"] = json.dumps(record.model_dump(mode="json"))
```

### Lineage Trace Walker
```python
# In sc_tools/provenance/trace.py
from pathlib import Path

def trace_lineage(file_path: str | Path, project_dir: str | Path = ".") -> list[dict]:
    """Walk provenance chain from file back to origins (D-11, D-12).

    Returns chronological list (oldest first) of lineage steps.
    """
    from sc_tools.provenance.sidecar import read_sidecar

    steps = []
    visited: set[str] = set()
    queue = [Path(file_path)]

    while queue:
        current = queue.pop(0)
        canonical = str(current.resolve())
        if canonical in visited:
            continue
        visited.add(canonical)

        prov = read_sidecar(current)
        # Fallback to adata.uns if sidecar missing and file is h5ad
        if prov is None and str(current).endswith(".h5ad"):
            prov = _read_uns_provenance(current)

        if prov is None:
            steps.append({
                "file": str(current),
                "command": None,
                "timestamp": None,
                "note": "origin (no provenance)",
            })
            continue

        steps.append({
            "file": str(current),
            "command": prov.get("command"),
            "timestamp": prov.get("timestamp"),
            "note": None,
        })

        for inp in prov.get("inputs", []):
            inp_path = Path(project_dir) / inp["path"]
            if not inp_path.exists():
                # D-10: SHA256 search for relocated files
                relocated = _find_by_sha256(inp["sha256"], project_dir)
                if relocated:
                    inp_path = relocated
                    steps[-1]["note"] = "input relocated"
            queue.append(inp_path)

    # D-12: chronological order (oldest first)
    steps.reverse()
    return steps
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| No provenance tracking | Minimal Provenance (Phase 2) | Phase 2 (2026-03) | command + timestamp + version only |
| Manual reproducibility notes | Automated sidecar writing | Phase 5 (this phase) | Full lineage DAG traversal |
| Hardcoded random seeds | Configurable random_state propagation | Phase 5 (this phase) | Deterministic clustering |

**Deprecated/outdated:**
- The minimal `Provenance` class stays in CLIResult for backwards compat but `ProvenanceRecord` is the full model for sidecars.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest 9.0.2 |
| Config file | pyproject.toml `[tool.pytest.ini_options]` |
| Quick run command | `pytest sc_tools/tests/test_provenance.py -x` |
| Full suite command | `pytest sc_tools/tests/ -x --tb=short` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| PRV-01 | Sidecar .provenance.json written alongside artifacts | unit | `pytest sc_tools/tests/test_provenance.py::test_sidecar_written_on_success -x` | Wave 0 |
| PRV-01 | No sidecar on error/no artifacts | unit | `pytest sc_tools/tests/test_provenance.py::test_no_sidecar_on_error -x` | Wave 0 |
| PRV-02 | Sidecar contains all required fields | unit | `pytest sc_tools/tests/test_provenance.py::test_sidecar_fields -x` | Wave 0 |
| PRV-02 | SHA256 checksum correctness | unit | `pytest sc_tools/tests/test_provenance.py::test_sha256_checksum -x` | Wave 0 |
| PRV-02 | Peak memory and runtime recorded | unit | `pytest sc_tools/tests/test_provenance.py::test_runtime_and_memory -x` | Wave 0 |
| PRV-03 | `sct provenance show` displays sidecar | integration | `pytest sc_tools/tests/test_cli_provenance.py::test_provenance_show -x` | Wave 0 |
| PRV-03 | `sct provenance show` falls back to adata.uns | integration | `pytest sc_tools/tests/test_cli_provenance.py::test_provenance_show_uns_fallback -x` | Wave 0 |
| PRV-04 | `sct provenance trace` walks lineage | integration | `pytest sc_tools/tests/test_cli_provenance.py::test_provenance_trace -x` | Wave 0 |
| PRV-04 | Trace handles missing intermediates gracefully | integration | `pytest sc_tools/tests/test_cli_provenance.py::test_provenance_trace_missing_input -x` | Wave 0 |
| PRV-05 | random_state threaded through _leiden_cluster | unit | `pytest sc_tools/tests/test_provenance.py::test_leiden_random_state -x` | Wave 0 |
| PRV-05 | Identical params produce identical clusters | unit | `pytest sc_tools/tests/test_provenance.py::test_leiden_reproducibility -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_provenance.py -x`
- **Per wave merge:** `pytest sc_tools/tests/ -x --tb=short`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/test_provenance.py` -- covers PRV-01, PRV-02, PRV-05 (unit tests for sidecar writing, checksums, Leiden reproducibility)
- [ ] `sc_tools/tests/test_cli_provenance.py` -- covers PRV-03, PRV-04 (CLI integration tests for show/trace commands)
- [ ] `sc_tools/tests/conftest.py` -- may need new fixtures: `adata_100_with_sidecar`, multi-step lineage fixture

## Open Questions

1. **How to track input files across all commands without modifying Phase 2-4 code?**
   - What we know: Commands already store `data["input"]` (single file) or load from `--from-dir`. The `artifacts` field captures outputs.
   - What is unclear: Whether to add `_input_files` convention key to data dict (minimal change) or add `input_files` field to CLIResult (cleaner but breaks compatibility).
   - Recommendation: Use `_input_files` convention key in data dict for now. Commands that produce artifacts already know their inputs. The cli_handler reads this key, builds InputFile records, and removes it before JSON emission. This avoids modifying CLIResult schema or touching Phase 2-4 code.

2. **adata.uns storage format for provenance**
   - What we know: h5ad stores `uns` values. Complex dicts need special handling.
   - What is unclear: Whether to store as JSON string (simple, portable) or as nested h5py groups (native but complex).
   - Recommendation: Store as JSON string (`json.dumps(record.model_dump(mode="json"))`). Simple to read back with h5py, no complex group structure, and round-trips cleanly.

3. **Depth limit for trace recursion**
   - What we know: Typical pipelines have 4-6 steps (raw -> QC -> preprocess -> integrate -> benchmark -> celltype). Cycle detection via visited set prevents infinite loops.
   - What is unclear: Whether very deep chains could cause issues.
   - Recommendation: No depth limit. Use visited-set cycle detection only. Real chains are shallow. Add a `--max-depth` option only if needed later.

## Sources

### Primary (HIGH confidence)
- `sc_tools/models/result.py` -- existing Provenance and CLIResult models (read directly)
- `sc_tools/cli/__init__.py` -- cli_handler decorator, _emit(), _state (read directly)
- `sc_tools/cli/discovery.py` -- register_discovery(app) pattern for command registration (read directly)
- `sc_tools/bm/integration.py` -- _leiden_cluster, compute_integration_metrics (read directly)
- `sc_tools/pp/reduce.py` -- leiden(), cluster() functions (read directly)
- Python stdlib docs: hashlib, resource, time.monotonic -- verified available on Python 3.10.15

### Secondary (MEDIUM confidence)
- `resource.getrusage` macOS vs Linux units -- verified empirically on target machine (macOS returns bytes, confirmed ru_maxrss = 10665984 for ~10MB process)
- scanpy.tl.leiden `random_state` parameter -- standard parameter in scanpy API, but cross-version determinism is best-effort

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all stdlib, no new dependencies
- Architecture: HIGH -- follows established cli_handler and register_discovery patterns from Phases 2-4
- Pitfalls: HIGH -- verified macOS/Linux ru_maxrss difference empirically, h5py access pattern well-understood
- Leiden reproducibility: MEDIUM -- random_state is necessary but may not guarantee bit-identical results across leidenalg versions

**Research date:** 2026-03-23
**Valid until:** 2026-04-23 (stable -- stdlib + existing patterns)
