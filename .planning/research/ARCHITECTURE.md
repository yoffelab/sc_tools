# Architecture Research

**Domain:** Agent-native CLI layer for single-cell/spatial transcriptomics library
**Researched:** 2026-03-20
**Confidence:** HIGH

## System Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                     CLI Surface (sct)                                │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐            │
│  │ sct qc   │  │ sct pp   │  │ sct bm   │  │ sct info │  ...       │
│  └────┬─────┘  └────┬─────┘  └────┬─────┘  └────┬─────┘            │
│       │              │             │              │                  │
├───────┴──────────────┴─────────────┴──────────────┴──────────────────┤
│                   Command Layer (thin adapters)                      │
│  ┌────────────────────────────────────────────────────────────────┐  │
│  │  CommandResult → JSON envelope {status, data, provenance}      │  │
│  └────────────────────────────────────────────────────────────────┘  │
│  ┌────────────┐  ┌─────────────┐  ┌────────────────┐                │
│  │ IO Gateway │  │  Provenance │  │ Output Formatter│               │
│  │ (load/save)│  │  Tracker    │  │ (JSON/human)   │                │
│  └─────┬──────┘  └──────┬──────┘  └───────┬────────┘                │
├────────┴─────────────────┴────────────────┴──────────────────────────┤
│                   sc_tools Library (unchanged)                       │
│  ┌─────────┐  ┌──────┐  ┌──────┐  ┌──────┐  ┌──────┐  ┌──────┐    │
│  │ ingest  │  │  qc  │  │  pp  │  │  tl  │  │  bm  │  │  pl  │    │
│  └─────────┘  └──────┘  └──────┘  └──────┘  └──────┘  └──────┘    │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐            │
│  │ storage  │  │ validate │  │ pipeline │  │ registry │            │
│  └──────────┘  └──────────┘  └──────────┘  └──────────┘            │
├─────────────────────────────────────────────────────────────────────┤
│                   Data Layer                                         │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐               │
│  │ .h5ad files  │  │ JSON sidecar │  │ registry.db  │               │
│  │ (AnnData)    │  │ (provenance) │  │ (optional)   │               │
│  └──────────────┘  └──────────────┘  └──────────────┘               │
└─────────────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

| Component | Responsibility | Typical Implementation |
|-----------|----------------|------------------------|
| CLI Surface (`sct`) | Parse args, validate inputs, dispatch to command layer | Typer app with subcommand groups (`sct qc`, `sct pp`, etc.) |
| Command Layer | Thin adapters: load data, call library, capture result, write provenance | One function per CLI command; returns `CommandResult` dataclass |
| IO Gateway | Memory-safe loading of h5ad files (backed mode, h5py selective reads) | Wraps `sc_tools.storage` with memory budget awareness |
| Provenance Tracker | Write JSON sidecar files alongside every output | Decorator/context manager that records inputs, params, outputs, timing |
| Output Formatter | Serialize `CommandResult` to JSON (default) or human-readable | `--human` flag switches from JSON to Rich tables/text |
| sc_tools Library | All actual computation (unchanged) | Existing modules: pp, qc, tl, bm, ingest, pl |
| Data Layer | Persistent storage of AnnData checkpoints, provenance, registry | h5ad files + `.provenance.json` sidecars |

## Recommended Project Structure

```
sc_tools/
├── cli/                        # NEW: CLI layer
│   ├── __init__.py             # Typer app factory, global options
│   ├── _types.py               # CommandResult, ProvenanceRecord dataclasses
│   ├── _output.py              # JSON/human output formatting
│   ├── _io.py                  # Memory-safe h5ad loading (backed mode, h5py)
│   ├── _provenance.py          # JSON sidecar writer, @track_provenance decorator
│   ├── qc.py                   # sct qc filter, sct qc report, sct qc metrics
│   ├── pp.py                   # sct pp preprocess, sct pp normalize, sct pp integrate
│   ├── ingest.py               # sct ingest load, sct ingest manifest
│   ├── bm.py                   # sct bm run, sct bm report, sct bm compare
│   ├── tl.py                   # sct tl score, sct tl annotate, sct tl deconvolve
│   ├── info.py                 # sct info, sct list-commands, sct describe <cmd>
│   └── assembly.py             # sct assembly build-mudata (late phase)
├── pp/                         # UNCHANGED
├── qc/                         # UNCHANGED
├── tl/                         # UNCHANGED
├── bm/                         # UNCHANGED
├── ingest/                     # UNCHANGED
├── pl/                         # UNCHANGED
├── storage.py                  # UNCHANGED
├── validate.py                 # UNCHANGED
├── pipeline.py                 # UNCHANGED
├── registry.py                 # UNCHANGED
└── ...
```

### Structure Rationale

- **`cli/` as separate package:** CLI is a pure interface layer. Isolating it prevents library code from depending on CLI concerns (arg parsing, output formatting). The library stays importable without Typer installed.
- **One file per domain:** Mirrors existing sc_tools module structure (qc, pp, tl, bm, ingest). Agents and humans navigate the same mental model.
- **Private `_types.py`, `_output.py`, `_io.py`, `_provenance.py`:** Cross-cutting concerns shared by all command files. Prefixed with `_` to signal internal-only.
- **No changes to existing library code:** The CLI wraps sc_tools; it does not modify it. This is a critical constraint from PROJECT.md.

## Architectural Patterns

### Pattern 1: Thin Adapter Commands

**What:** Each CLI command is a thin function that (1) loads data via IO Gateway, (2) calls exactly one sc_tools function, (3) wraps the result in CommandResult, (4) writes provenance sidecar. No business logic in the CLI layer.

**When to use:** Every command. This is the universal pattern.

**Trade-offs:** Pro: library and CLI evolve independently, library stays testable without CLI. Con: some boilerplate per command. Worth it because the boilerplate is mechanical and keeps boundaries clean.

**Example:**
```python
# sc_tools/cli/qc.py
import typer
from sc_tools.cli._types import CommandResult
from sc_tools.cli._io import load_adata_safe
from sc_tools.cli._provenance import track_provenance
from sc_tools.cli._output import emit

app = typer.Typer(name="qc", help="Quality control operations")

@app.command()
@track_provenance
def filter(
    input: str = typer.Argument(..., help="Path to h5ad file"),
    output: str = typer.Option("filtered.h5ad", help="Output path"),
    min_genes: int = typer.Option(200),
    min_counts: int = typer.Option(500),
    human: bool = typer.Option(False, help="Human-readable output"),
) -> None:
    """Filter low-quality cells from AnnData."""
    from sc_tools.qc import calculate_qc_metrics, filter_cells

    adata = load_adata_safe(input)
    calculate_qc_metrics(adata)
    n_before = adata.n_obs
    filter_cells(adata, min_genes=min_genes, min_counts=min_counts)
    adata.write_h5ad(output)

    result = CommandResult(
        status="success",
        data={"n_before": n_before, "n_after": adata.n_obs, "output": output},
    )
    emit(result, human=human)
```

### Pattern 2: CommandResult Envelope

**What:** Every command returns the same JSON envelope structure. Agents parse `status` and `data` fields; humans get formatted tables. Errors follow the same envelope with `status: "error"`.

**When to use:** All commands. Consistency is the point.

**Trade-offs:** Pro: agents always know the response shape; error handling is uniform. Con: slightly verbose for trivial commands. Worth it because agents parse this programmatically.

**Example:**
```python
# sc_tools/cli/_types.py
from dataclasses import dataclass, field, asdict
from typing import Any
import json

@dataclass
class CommandResult:
    status: str  # "success" | "error"
    data: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] | None = None
    warnings: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, default=str)
```

**JSON output (stdout):**
```json
{
  "status": "success",
  "data": {
    "n_before": 125000,
    "n_after": 98432,
    "output": "results/adata.filtered.h5ad"
  },
  "provenance": {
    "command": "sct qc filter",
    "input": "results/adata.raw.h5ad",
    "output": "results/adata.filtered.h5ad",
    "params": {"min_genes": 200, "min_counts": 500},
    "duration_seconds": 12.3,
    "timestamp": "2026-03-20T14:30:00Z"
  },
  "warnings": []
}
```

### Pattern 3: JSON Sidecar Provenance

**What:** Every output file gets a companion `.provenance.json` sidecar. The sidecar records the command, inputs (with checksums), parameters, outputs, timing, and environment. This is file-based lineage -- no database required.

**When to use:** Every command that produces an output file.

**Trade-offs:** Pro: zero infrastructure, grep-able, version-controllable, schema can evolve freely. Con: not queryable like a DB (addressed later by optional DB import). This matches the PROJECT.md decision: "file-based provenance before DB."

**Example sidecar (`results/adata.filtered.h5ad.provenance.json`):**
```json
{
  "command": "sct qc filter",
  "version": "0.4.0",
  "timestamp": "2026-03-20T14:30:00Z",
  "duration_seconds": 12.3,
  "inputs": [
    {"path": "results/adata.raw.h5ad", "sha256": "abc123..."}
  ],
  "outputs": [
    {"path": "results/adata.filtered.h5ad", "sha256": "def456..."}
  ],
  "params": {"min_genes": 200, "min_counts": 500},
  "environment": {
    "hostname": "brb-gpu01",
    "gpu": "A100-80G",
    "peak_memory_gb": 18.2
  }
}
```

### Pattern 4: Memory-Safe IO Gateway

**What:** A wrapper around `sc_tools.storage` that chooses the loading strategy based on file size and available memory. For files under a threshold (configurable, default 8GB), load fully. For larger files, use h5py selective reads or AnnData backed mode. For metrics-only operations, subsample.

**When to use:** Every command that reads h5ad files.

**Trade-offs:** Pro: prevents OOM on 25G Visium HD files. Con: backed mode limits which operations are possible (only X is mutable). Mitigated by detecting which operations need full load vs. selective access.

**Example:**
```python
# sc_tools/cli/_io.py
import os
from pathlib import Path

def load_adata_safe(path: str, max_memory_gb: float = 8.0, subsample: int | None = None):
    """Load AnnData with memory-aware strategy selection."""
    import anndata as ad

    file_size_gb = Path(path).stat().st_size / (1024 ** 3)

    if subsample is not None:
        # For metrics: load obs, subsample indices, load X[indices]
        return _load_subsampled(path, n=subsample)
    elif file_size_gb <= max_memory_gb:
        return ad.read_h5ad(path)
    else:
        # Backed mode for large files
        return ad.read_h5ad(path, backed="r")

def _load_subsampled(path: str, n: int):
    """Load only n random cells via h5py for metrics computation."""
    import h5py
    import numpy as np
    import anndata as ad

    with h5py.File(path, "r") as f:
        n_obs = f["X"].shape[0]
        indices = np.sort(np.random.choice(n_obs, min(n, n_obs), replace=False))
        # Build AnnData from selected rows only
        ...
    return adata
```

### Pattern 5: Self-Describing Discovery

**What:** `sct list-commands` returns a JSON array of all available commands with their descriptions, required inputs, and output types. `sct describe <command>` returns the full schema for a command. This lets agents discover capabilities at runtime.

**When to use:** Agent integration. Agents call `sct list-commands` to know what is available, then `sct describe qc.filter` to get the parameter schema.

**Trade-offs:** Pro: agents never need hardcoded knowledge of the CLI surface. Con: requires maintaining schema metadata alongside commands. Mitigated by deriving schemas from Typer's type annotations.

## Data Flow

### Command Execution Flow

```
Agent/User invokes: sct pp preprocess --input adata.raw.h5ad --modality visium
    |
    v
[Typer CLI] parse args, validate types
    |
    v
[IO Gateway] check file size (25G) -> backed mode or subsample
    |
    v
[Provenance Tracker] start timer, record inputs + params
    |
    v
[sc_tools.pp.preprocess()] actual computation (unchanged library code)
    |
    v
[IO Gateway] write output h5ad
    |
    v
[Provenance Tracker] stop timer, compute output checksum, write sidecar JSON
    |
    v
[Output Formatter] serialize CommandResult to JSON on stdout
    |
    v
Agent reads JSON from stdout -> decides next step
```

### Provenance Chain

```
adata.raw.h5ad                    # Phase 0 output
    |
    +-- adata.raw.h5ad.provenance.json
    |
    v
sct qc filter --input adata.raw.h5ad --output adata.filtered.h5ad
    |
    v
adata.filtered.h5ad               # Phase 1 output
    |
    +-- adata.filtered.h5ad.provenance.json  (links to adata.raw.h5ad via inputs)
    |
    v
sct pp preprocess --input adata.filtered.h5ad --output adata.integrated.h5ad
    |
    v
adata.integrated.h5ad             # Phase 3 output
    |
    +-- adata.integrated.h5ad.provenance.json  (links to adata.filtered.h5ad)
```

Each sidecar's `inputs[].path` field forms a DAG that can be walked to reconstruct the full lineage without a database.

### Key Data Flows

1. **Ingestion flow:** Raw vendor output -> `sct ingest load` -> h5ad + sidecar. Multiple samples -> `sct ingest concat` -> merged h5ad + sidecar listing all input sidecars.
2. **Analysis flow:** h5ad -> `sct qc filter` -> `sct pp preprocess` -> `sct tl annotate` -> final h5ad. Each step reads previous output, writes new output + sidecar.
3. **Benchmarking flow:** Multiple h5ad files (one per method) -> `sct bm compare` -> metrics JSON + report HTML. Reads selectively via h5py (no full load).
4. **Discovery flow:** Agent calls `sct list-commands` -> gets JSON array -> selects command -> calls `sct describe <cmd>` -> gets param schema -> constructs and executes command.

## Scaling Considerations

| Scale | Architecture Adjustments |
|-------|--------------------------|
| Small datasets (<500K cells, <4G) | Default full-memory load. No special handling needed. |
| Medium datasets (500K-2M cells, 4-15G) | IO Gateway auto-selects backed mode. Most operations work. Some methods (scVI) need full load -- subsample first. |
| Large datasets (2M+ cells, 15-25G) | Mandatory backed mode or h5py selective reads. Metrics via subsampling (50K-100K cells). Integration on HPC with GPU. CLI generates SLURM scripts via `sct slurm generate`. |

### Scaling Priorities

1. **First bottleneck: Memory on load.** The 25G Visium HD files cannot be fully loaded on most machines. The IO Gateway's memory-budget strategy is the critical feature. Build this in Phase 1.
2. **Second bottleneck: Integration runtime.** scVI on 2.5M cells takes hours even on GPU. The CLI does not solve this directly -- it delegates to HPC via SLURM script generation (`sct slurm generate`), which already exists in the library.
3. **Third bottleneck: Benchmarking I/O.** Loading multiple 25G h5ad files for comparison. Solved by h5py selective reads -- load only `obsm['X_scvi']` and `obs` columns needed for metrics, never full X.

## Anti-Patterns

### Anti-Pattern 1: Business Logic in CLI Layer

**What people do:** Put data transformations, filtering logic, or algorithmic decisions in CLI command functions.
**Why it's wrong:** CLI commands become untestable without CLI invocation. Library and CLI couple tightly. Breaks the "CLI wraps sc_tools, not replaces it" constraint.
**Do this instead:** CLI commands call exactly one library function. If you need a new operation, add it to the library first, then wrap it.

### Anti-Pattern 2: Eager Full Load of Large Files

**What people do:** `adata = sc.read_h5ad(path)` for every command, even when only obs metadata is needed.
**Why it's wrong:** OOM on Visium HD datasets. The PROJECT.md origin story is literally about an agent script that OOM'd at 44G because it loaded everything.
**Do this instead:** IO Gateway checks file size, selects strategy (full/backed/h5py-selective/subsample). Commands declare their access pattern (full, obs-only, obsm-selective).

### Anti-Pattern 3: Unstructured stdout Mixing

**What people do:** `print()` status messages, progress, and results all to stdout.
**Why it's wrong:** Agents parse stdout as JSON. Any non-JSON text breaks parsing. Progress messages corrupt the output stream.
**Do this instead:** JSON result on stdout. All logging, progress, and status messages go to stderr. Use `typer.echo(msg, err=True)` for human-visible messages.

### Anti-Pattern 4: Provenance as Afterthought

**What people do:** Build all commands first, add provenance tracking later as a separate pass.
**Why it's wrong:** Retrofitting provenance means inconsistent coverage, missed edge cases, and schema that does not match actual usage patterns.
**Do this instead:** Build provenance tracking in Phase 1 as a decorator/context manager. Every command gets it from day one via `@track_provenance`.

### Anti-Pattern 5: Fat Command Objects

**What people do:** Create elaborate Command/Handler/Executor class hierarchies with abstract base classes and factory patterns.
**Why it's wrong:** Over-engineering for a CLI that wraps an existing library. Each command is fundamentally: load, call, save, report.
**Do this instead:** Plain functions with a Typer decorator. The `CommandResult` dataclass is the only shared abstraction. Keep it flat.

## Integration Points

### External Services

| Service | Integration Pattern | Notes |
|---------|---------------------|-------|
| HPC (SLURM) | `sct slurm generate` produces sbatch scripts | Delegates to existing `sc_tools/ingest/slurm.py` and `sc_tools/bm/slurm.py` |
| Cloud storage (S3/GCS) | IO Gateway uses `sc_tools.storage.resolve_fs()` | fsspec handles all URI resolution; CLI just passes paths through |
| Registry DB (SQLAlchemy) | Optional: `sct registry status` reads from DB | Existing `sc_tools/registry.py`; CLI wraps `_cli_status()` |
| MCP server | Coexists; MCP calls library directly, CLI is for subprocess invocation | Agents may use MCP for in-process calls or CLI for subprocess isolation |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| CLI -> Library | Direct Python function calls | CLI imports and calls library functions. No serialization boundary. |
| CLI -> IO Gateway | Function calls with path + memory budget | IO Gateway wraps `sc_tools.storage` with memory-awareness |
| CLI -> Provenance | Decorator + context manager | `@track_provenance` wraps command function; writes sidecar on completion |
| CLI -> Output Formatter | `emit(result, human=bool)` | Single function that serializes to JSON (stdout) or Rich table (stdout) |
| Command files -> _types | Import `CommandResult` | Shared dataclass; all commands return the same envelope |

### CLI vs MCP: When to Use Each

| Scenario | Use CLI | Use MCP |
|----------|---------|---------|
| Agent runs analysis step | Yes -- subprocess isolation, structured JSON output | Also works -- in-process, lower overhead |
| Snakemake rule calls sc_tools | Yes -- Snakemake invokes CLI as shell command | No -- MCP requires running server |
| Human runs analysis | Yes -- with `--human` flag | No -- MCP is agent-only |
| Claude Code orchestrator | Either -- MCP for quick queries, CLI for heavy compute | Either |
| HPC batch job | Yes -- CLI is the entry point in sbatch scripts | No -- MCP server not running on compute nodes |

## Build Order (Dependency Graph)

The CLI layer should be built in this order, where each phase depends on the previous:

```
Phase 1: Foundation
    _types.py (CommandResult)
    _output.py (JSON/human formatter)
    _io.py (memory-safe loading)
    _provenance.py (sidecar writer + decorator)
    __init__.py (Typer app factory, global --human/--verbose flags)
    info.py (sct list-commands, sct describe)
        |
        v
Phase 2: Core Commands (highest immediate value)
    qc.py (filter, metrics, report -- most common operations)
    pp.py (preprocess, normalize, integrate)
    ingest.py (load, manifest, concat)
        |
        v
Phase 3: Analysis + Benchmarking
    tl.py (score, annotate, deconvolve)
    bm.py (run, compare, report -- replaces existing bm/cli.py)
        |
        v
Phase 4: Multi-modal Assembly
    assembly.py (build-mudata, cross-modal queries)
```

**Rationale:** Phase 1 infrastructure (types, IO, provenance) must exist before any command. Phase 2 covers the commands agents call most (QC and preprocessing). Phase 3 adds analysis. Phase 4 is late-stage because MuData assembly requires all modalities processed independently first.

## Technology Decision: Typer over argparse

**Use Typer** (built on Click) for the CLI framework because:

1. **Type-hint driven:** Function signatures become the CLI schema. No separate argparse definitions to maintain. This is critical for self-describing discovery -- `sct describe` can introspect Typer commands to generate parameter schemas.
2. **Subcommand composition:** `app.add_typer(qc_app, name="qc")` mirrors the sc_tools module structure naturally.
3. **Testing:** `typer.testing.CliRunner` captures output without subprocess overhead.
4. **Rich integration:** Human-readable output gets tables, colors, progress bars for free via `--human` flag.
5. **Click compatibility:** Existing `bm/cli.py` uses argparse; migration to Typer is straightforward since Typer wraps Click which is argparse-compatible in spirit.

The existing `bm/cli.py` (argparse-based, 263 lines) demonstrates the pain: manual parser construction, manual dispatch, no structured output, no provenance. Typer eliminates the boilerplate while adding type safety.

## Sources

- [Nextflow CLI architecture and command dispatch](https://deepwiki.com/nextflow-io/nextflow/2.2-command-line-interface)
- [Nextflow workflow outputs: JSON params in, JSON manifest out](https://seqera.io/podcasts/episode-57-pipeline-chaining-meta-pipelines-part-2/)
- [Snakemake codebase architecture: CLI -> API -> engine](https://snakemake.readthedocs.io/en/stable/project_info/codebase.html)
- [Cell Ranger CLI: verb-based subcommands with structured outs/ directory](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments)
- [Typer: type-hint-driven CLI framework](https://github.com/fastapi/typer)
- [scanpy.read_h5ad backed mode for memory efficiency](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html)
- [scanpy backed mode limitations with large datasets](https://github.com/scverse/scanpy/issues/2365)
- [HyProv: sidecar-based provenance for scientific workflows](https://arxiv.org/html/2511.07574)
- [yProv4ML: PROV-JSON format for ML provenance](https://arxiv.org/pdf/2507.01075)
- [Python CLI framework comparison 2025: Click at 38.7% adoption](https://dasroot.net/posts/2025/12/building-cli-tools-python-click-typer-argparse/)

---
*Architecture research for: agent-native CLI layer on sc_tools*
*Researched: 2026-03-20*
