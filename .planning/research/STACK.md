# Technology Stack: Agent-Native CLI Layer

**Project:** sc_tools Agent-Native CLI
**Researched:** 2026-03-20
**Overall confidence:** HIGH

## Recommended Stack

### CLI Framework

| Technology | Version | Purpose | Why | Confidence |
|------------|---------|---------|-----|------------|
| Typer | >= 0.24.1 | CLI framework | Type-hint-driven, auto-generates help, built on Click, Rich integration for human output. Dominant modern Python CLI framework (38.7% of new projects use Click; Typer inherits that ecosystem while adding type safety). Existing sc_tools uses Python 3.11+ type hints everywhere -- Typer leverages those directly as CLI parameters. | HIGH |
| Rich | >= 13.0 | Human-readable output | Already a Typer dependency. Use `rich.console.Console(stderr=True)` for human output, keeping stdout clean for JSON. Tables, progress bars, syntax highlighting for `--human` mode. | HIGH |

**Why not Click directly:** Click requires explicit decorator-based parameter definitions. Typer auto-derives parameters from function signatures with type hints -- less boilerplate, and the function signatures double as the self-describing schema (see below). sc_tools functions already have typed signatures.

**Why not argparse:** No auto-completion, no rich help formatting, verbose boilerplate. The existing scripts use argparse (e.g., `run_qc_report.py`, `run_preprocessing.py`) -- migrating to Typer reduces LOC and adds shell completion for free.

**Why not Fire:** Google Fire is convenient for quick prototyping but has poor control over help text, no shell completion, and unreliable type coercion. Not suitable for agent-facing tools that need predictable structured output.

### Structured Output

| Technology | Version | Purpose | Why | Confidence |
|------------|---------|---------|-----|------------|
| Pydantic | >= 2.7 | Output schemas | Define result models as Pydantic BaseModel subclasses. `.model_dump_json()` for agent output, `.model_dump()` for internal use. Typer already depends on Pydantic (via Click parameter types). Pydantic v2 is fast (Rust core). Already in the Python ecosystem for sc_tools' MCP server. | HIGH |
| stdlib json | 3.11+ | JSON serialization | For simple cases where Pydantic is overkill. `json.dumps()` to stdout. | HIGH |

**Output format pattern:** JSON to stdout (default, for agents), human-readable via Rich to stderr when `--human` flag is set. This follows the Unix convention that tools like `jq`, `rg --json`, `gh`, and AWS CLI all use.

**Why not JSONL:** JSONL (newline-delimited JSON) is useful for streaming/log output, but most sc_tools operations are batch (run, return result). Use JSONL only for progress events on long-running operations (e.g., integration benchmarks). Standard JSON for command results.

**Why not YAML:** More human-readable but harder to parse reliably, slower to serialize, and not what agents expect. JSON is the universal agent interchange format.

**Why not Protocol Buffers / MessagePack:** Overkill for a CLI tool. Binary formats don't pipe well. JSON is universal and debuggable.

### Output Schema Contract

```python
# Every CLI command returns one of these to stdout:
from pydantic import BaseModel
from typing import Any

class CLIResult(BaseModel):
    """Standard envelope for all sct command output."""
    status: str  # "ok" | "error" | "warning"
    command: str  # e.g., "sct qc report"
    data: dict[str, Any]  # command-specific payload
    artifacts: list[str]  # file paths created/modified
    provenance: dict[str, Any] | None = None  # optional sidecar ref
    message: str = ""  # human-readable summary
```

This envelope pattern lets agents reliably parse any `sct` command output: check `status`, extract `data`, find `artifacts`.

### Provenance Tracking

| Technology | Version | Purpose | Why | Confidence |
|------------|---------|---------|-----|------------|
| Custom JSON sidecars | n/a | Lightweight provenance | Write a `.provenance.json` file alongside each output artifact. Simple, no dependencies, schema evolves with usage. Matches PROJECT.md decision: "file-based first, DB later." | HIGH |
| Pydantic | >= 2.7 | Sidecar schema validation | Define provenance record as a Pydantic model. Validates on write, serializes cleanly. | HIGH |
| prov (W3C PROV) | >= 2.1.1 | **Not recommended yet** | Full W3C PROV-DM implementation. Powerful but heavyweight for the "let the schema emerge" phase. Revisit when provenance model stabilizes and needs interoperability with external systems. | MEDIUM |
| ro-crate-py | >= 0.12 | **Not recommended yet** | RO-Crate is the gold standard for FAIR research packaging, used by Galaxy, WorkflowHub, CWL. But it solves a different problem (packaging for publication/sharing), not runtime lineage tracking. Consider for export/archival phase later. | MEDIUM |
| provit | ~0.3 | **Not recommended** | JSON-LD sidecar approach (inspiration for our design), but unmaintained (last release ~2021), limited adoption. Take the idea, not the dependency. | LOW |

**Provenance sidecar format:**

```json
{
  "schema_version": "0.1.0",
  "created": "2026-03-20T14:30:00Z",
  "command": "sct preprocess --method harmony --adata input.h5ad",
  "inputs": [
    {"path": "results/adata.filtered.h5ad", "checksum": "sha256:abc123..."}
  ],
  "outputs": [
    {"path": "results/preprocessing/harmony.h5ad", "checksum": "sha256:def456..."}
  ],
  "parameters": {
    "method": "harmony",
    "n_hvgs": 3000,
    "batch_key": "sample"
  },
  "environment": {
    "sc_tools_version": "0.1.0",
    "python": "3.11.8",
    "gpu": false
  },
  "duration_seconds": 142.5,
  "agent": "claude-code"
}
```

**Why this over W3C PROV / RO-Crate:** The project decision is explicit: "file-based provenance before DB." The schema needs to emerge from actual usage patterns. Starting with a custom Pydantic model means: (1) zero new dependencies, (2) schema changes are just model changes, (3) easy to migrate to W3C PROV or RO-Crate later by writing an export function. Starting with a formal ontology before the domain model stabilizes is premature abstraction.

**Migration path:** Custom JSON sidecars -> W3C PROV-JSON export (add `prov` library) -> RO-Crate packaging for publication (add `rocrate`). Each step is additive, not a rewrite.

### Self-Describing CLI

| Technology | Version | Purpose | Why | Confidence |
|------------|---------|---------|-----|------------|
| Typer introspection | >= 0.24.1 | Command discovery | Typer exposes its command tree via Click's internal API. Build `sct list-commands` and `sct describe <cmd>` by walking the Click Group tree at runtime. | HIGH |
| Custom JSON schema export | n/a | Machine-readable command catalog | A `sct schema` command that exports all commands, their parameters, types, defaults, and descriptions as JSON. Agents call this once to discover what's available. | HIGH |
| Pydantic | >= 2.7 | Schema generation | `CLIResult.model_json_schema()` exports the output schema. Combined with parameter introspection, gives agents a full contract. | HIGH |

**Why not OpenAPI:** OpenAPI is for HTTP APIs, not CLIs. There is no widely-adopted equivalent standard for CLI tools. The closest patterns are: (1) `gh` (GitHub CLI) uses `--json` flag + field selection, (2) AWS CLI uses `--output json` + JMESPath queries, (3) `kubectl` uses `-o json` + JSONPath. We follow this precedent with a simpler approach: export a JSON schema that agents parse.

**Self-description pattern:**

```bash
# Agent discovers commands
sct schema                    # Full JSON schema of all commands + parameters + output types
sct list-commands             # Simple list: ["qc", "preprocess", "ingest", ...]
sct describe preprocess       # Detailed JSON: parameters, types, defaults, description

# Human discovers commands
sct --help                    # Typer auto-generated help (Rich-formatted)
sct preprocess --help         # Per-command help
```

The `sct schema` output would be a JSON document:

```json
{
  "name": "sct",
  "version": "0.1.0",
  "commands": {
    "preprocess": {
      "description": "Run preprocessing pipeline on AnnData",
      "parameters": {
        "adata": {"type": "path", "required": true, "description": "Input h5ad file"},
        "method": {"type": "string", "enum": ["harmony", "scvi", "resolvi"], "default": "harmony"},
        "n_hvgs": {"type": "integer", "default": 3000}
      },
      "output_schema": {"$ref": "#/definitions/PreprocessResult"}
    }
  }
}
```

### Supporting Libraries

| Library | Version | Purpose | When to Use | Confidence |
|---------|---------|---------|-------------|------------|
| pydantic | >= 2.7 | Schema validation, JSON serialization | Every command (output models, provenance records) | HIGH |
| typer | >= 0.24.1 | CLI framework | Entry point (`sct` command) | HIGH |
| rich | >= 13.0 | Human-readable output | When `--human` flag or interactive terminal | HIGH |
| hashlib (stdlib) | 3.11+ | File checksums for provenance | Provenance sidecar generation | HIGH |
| datetime (stdlib) | 3.11+ | Timestamps | Provenance records | HIGH |

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| CLI framework | Typer 0.24 | Click 8.1 | More boilerplate; Typer wraps Click anyway |
| CLI framework | Typer 0.24 | argparse (stdlib) | No completion, no rich help, verbose. Current scripts use it -- migration reduces LOC |
| CLI framework | Typer 0.24 | Fire 0.7 | Poor help text control, unreliable type coercion |
| Output format | JSON (Pydantic) | YAML | Slower, ambiguous parsing, not agent-native |
| Output format | JSON (Pydantic) | TOML | Not suitable for structured data output (config format) |
| Provenance | Custom JSON sidecars | W3C PROV (prov lib) | Premature formalization; add later when model stabilizes |
| Provenance | Custom JSON sidecars | RO-Crate | Solves packaging/sharing, not runtime tracking; add for export later |
| Provenance | Custom JSON sidecars | SQLite/registry DB | PROJECT.md explicitly defers DB provenance until model stabilizes |
| Self-description | JSON schema export | OpenAPI | Designed for HTTP APIs, not CLIs |
| Self-description | JSON schema export | man pages | Not machine-readable |

## Installation

```bash
# Add to pyproject.toml [project.optional-dependencies]
# cli extra (new)
pip install -e ".[cli]"

# Which installs:
# typer[all] >= 0.24.1   (includes rich, shellingham)
# pydantic >= 2.7         (already in ecosystem via mcp/scvi-tools)
```

The `[cli]` extra keeps the CLI layer optional -- existing library users don't need typer. The `pydantic` dependency is already satisfied transitively by `mcp`, `scvi-tools`, and other existing extras.

**pyproject.toml addition:**

```toml
[project.optional-dependencies]
cli = [
    "typer[all] >= 0.24.1",
    "pydantic >= 2.7",
]

[project.scripts]
sct = "sc_tools.cli:app"
```

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| JSON to stdout, human to stderr | Unix convention. Agents parse stdout; humans see Rich on stderr. No output corruption. |
| Pydantic envelope for all output | Agents always know: check `status`, read `data`, find `artifacts`. Predictable parsing. |
| `--human` flag (opt-in) not `--json` (opt-out) | Agent-native means JSON is default. Humans opt into readable mode. Inverts the traditional pattern. |
| Custom provenance before standards | Let domain model emerge from usage. W3C PROV / RO-Crate are migration targets, not starting points. |
| `sct schema` for self-description | One command gives agents the full contract. No need to parse help text. |
| Typer over Click | Same ecosystem (Typer wraps Click), but type-hint derivation matches sc_tools' existing code style. |
| Separate `[cli]` extra | CLI is an interface layer -- library users shouldn't need typer. Keeps core lean. |

## Version Compatibility

| Component | Min Version | Max Tested | Notes |
|-----------|-------------|------------|-------|
| Python | 3.11 | 3.14 | Matches existing sc_tools requirement |
| Typer | 0.24.1 | 0.24.1 | Current release as of 2026-02 |
| Pydantic | 2.7 | 2.x | v2 required (Rust core, `model_dump_json`) |
| Rich | 13.0 | latest | Typer[all] pins compatible version |

## Sources

- [Typer on PyPI](https://pypi.org/project/typer/) - v0.24.1, Python 3.10+, Feb 2026 (HIGH confidence)
- [Typer GitHub](https://github.com/fastapi/typer) - maintained by FastAPI team (HIGH confidence)
- [prov on PyPI](https://pypi.org/project/prov/) - v2.1.1, W3C PROV-DM implementation (HIGH confidence)
- [provit on GitHub](https://github.com/diggr/provit) - JSON-LD sidecar inspiration, unmaintained (LOW confidence)
- [ro-crate-py on GitHub](https://github.com/ResearchObject/ro-crate-py) - RO-Crate 1.2, Python 3.9+ (MEDIUM confidence)
- [RO-Crate PLOS ONE paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0309210) - Workflow Run RO-Crate profiles (MEDIUM confidence)
- [jc CLI tool](https://kellyjonbrazil.github.io/jc/) - JSON conversion patterns for CLI tools (HIGH confidence)
- [CLI framework comparison (2025)](https://dasroot.net/posts/2025/12/building-cli-tools-python-click-typer-argparse/) - Click 38.7% market share (MEDIUM confidence)
- [AWS CLI output formats](https://docs.aws.amazon.com/cli/latest/userguide/cli-usage-output-format.html) - `--output json` pattern precedent (HIGH confidence)
- [Pydantic AI agent patterns](https://ai.pydantic.dev/) - Pydantic for structured agent output (MEDIUM confidence)

---
*Stack analysis: 2026-03-20*

---

## v2.0 Milestone Addendum: Spatial Plots, UMAP Plots, HTML Reports, sct concat

**Scope:** NEW capabilities only for v2.0 (Report Plots & Sample Concat milestone).
**Researched:** 2026-03-27
**Confidence:** HIGH — all findings verified against existing codebase.

### Finding: No New Dependencies Required

Every library needed for v2.0 is already declared in `pyproject.toml`. The work is entirely
wiring and implementation, not procurement.

**One promotion needed:** `plotly` moves from `[pipeline]` optional to base `dependencies`.
**One maintenance fix:** hardcoded Plotly CDN version tag in `report_utils.py` is stale.

---

### Spatial Plots (log1p_counts, log1p_genes, %mt, pass/fail overlay)

**Library:** `scanpy` (base dep, `>=1.9`) — `sc.pl.spatial`
**Confidence:** HIGH — pattern already implemented in `sc_tools/qc/report.py` lines 314–333.

The post-filter report already generates `log1p_total_counts` spatial plots per sample using:

```python
sc.pl.spatial(sub, color=col, library_id=lib, show=False, ax=ax, frameon=False)
```

The v2.0 QC report extends this to 4 panels per sample in a `plt.subplots(1, 4)` grid:

| Panel | obs column | Derivation |
|-------|-----------|------------|
| log1p_counts | `log1p_total_counts` | `np.log1p(obs["total_counts"])` — already computed |
| log1p_genes | `log1p_n_genes_by_counts` | `np.log1p(obs["n_genes_by_counts"])` — add same way |
| %mt | `pct_counts_mt` | Already in obs when mt genes present — guard with `has_mt` |
| pass/fail | `qc_pass_str` | Cast `classified["qc_pass"]` bool to `"Pass"`/`"Fail"` string for color |

Celltype report spatial plots (per cell type per sample) use the same pattern with `color=celltype_col`.

**No new library. No version change needed.**

---

### UMAP Plots (leiden, batch, patient, bio obs columns)

**Library:** `scanpy` (base dep, `>=1.9`) — `sc.pl.umap`
**Confidence:** HIGH — pattern already used in `sc_tools/pl/benchmarking.py` line 529.

```python
sc.pl.umap(adata, color=col, show=False, ax=ax, frameon=False)
```

The integration report generates a UMAP grid over obs columns: `leiden_key`, `batch_key`,
`patient_col`, and any bio obs columns used in scib metrics. Caller ensures `X_umap` is present
(data precondition, not a library gap — the preprocess phase always computes UMAP).

**No new library. No version change needed.**

---

### HTML Embedding: Static Plots (spatial, UMAP)

**Mechanism:** `fig_to_base64(fig)` in `sc_tools/qc/report_utils.py` — converts matplotlib Figure
to base64-encoded PNG string embedded as `<img src="data:image/png;base64,...">`.

Static PNG is the correct choice for spatial and UMAP plots. Reasons:

1. **Scale**: Visium HD datasets reach 2.5M spots. An interactive Plotly scatter with 2.5M points
   makes the browser unusable. Static PNG at 150 dpi is sufficient for QC review.
2. **Self-contained**: Base64 PNG requires no CDN, no JS, no internet connection.
3. **Consistency**: All existing QC plots already use `fig_to_base64`. New plots should follow
   the same pattern for template consistency.

**No new library. No change needed.**

---

### HTML Embedding: Interactive Plots (Plotly)

**Mechanism:** `plotly_to_html(fig)` calls `fig.to_html(full_html=False, include_plotlyjs=False)`.
One CDN `<script>` tag is injected per report in `_wrap_with_tabs`.

**Action required — promote plotly to base deps:**

`plotly` is currently only in `[pipeline]` optional. The QC/integration/celltype report functions
already call `plotly_to_html` (radar charts, batch-vs-bio scatter in `compute_integration_section`).
Any user running `sct qc report` or `sct integration report` without `[pipeline]` extras will get
an ImportError at runtime. Promote to base:

```toml
# pyproject.toml [project.dependencies]
"plotly>=5.18",   # add to base dependencies (currently only in [pipeline])
```

`plotly>=5.18` supports Python 3.8+, so 3.10 compatibility is confirmed. The `[pipeline]`
entry can remain as a no-op duplicate (pip deduplicates).

**Action required — update stale CDN version tag:**

`report_utils.py` line 517 hardcodes `plotly-2.27.0.min.js` (a 2023 release). Current plotly.py
is 6.x. Use `plotly-latest.min.js` to track current stable automatically:

```python
# report_utils.py line 517 — change from:
plotly_cdn = '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
# to:
plotly_cdn = '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'
```

**Known limitation:** CDN-based reports require internet access. This is an existing design choice.
Offline mode (embedding the ~3MB plotly.js bundle via `include_plotlyjs=True`) is not worth the
file size for HPC/local use. Document as limitation in report footer.

---

### sct concat CLI Command

**Backend:** `sc_tools.ingest.concat_samples` — fully implemented in
`sc_tools/ingest/loaders.py` lines 903–947.

Signature:
```python
def concat_samples(
    adatas: list[ad.AnnData],
    *,
    sample_col: str = "sample",
    calculate_qc: bool = True,
) -> ad.AnnData
```

Uses `ad.concat(adatas, merge="same", uns_merge="same")` which preserves spatial coordinates.
Runs `sc.pp.calculate_qc_metrics` post-concat if `calculate_qc=True`.

**What v2.0 adds:** A CLI wrapper only. No changes to backend logic.

```python
# sc_tools/cli/commands/concat.py — sketch
@app.command("concat")
def concat_cmd(
    input_files: Annotated[list[Path], typer.Argument(...)],
    output: Annotated[Path, typer.Option(...)],
    sample_col: str = "library_id",
    calculate_qc: bool = True,
):
    from sc_tools.ingest import concat_samples
    import anndata as ad
    adatas = [ad.read_h5ad(f) for f in input_files]
    result = concat_samples(adatas, sample_col=sample_col, calculate_qc=calculate_qc)
    result.write_h5ad(output)
    # emit CLIResult JSON to stdout
```

For large files, use backed mode before concat to avoid OOM: `ad.read_h5ad(f, backed="r")` then
`.to_memory()` per-sample before passing to `concat_samples`. This is consistent with the IO
Gateway pattern from Phase 7 (memory safety).

**No new library. Uses existing anndata, scanpy, typer, pydantic.**

---

### v2.0 Stack Changes Summary

| Item | Action | pyproject.toml change |
|------|--------|-----------------------|
| `plotly` | Promote from `[pipeline]` optional to base `dependencies` | Add `"plotly>=5.18"` to `[project.dependencies]` |
| `scanpy` | No change — `sc.pl.spatial` + `sc.pl.umap` already available | None |
| `squidpy` | No change — already used, not needed for plot generation | None |
| `matplotlib` | No change — `fig_to_base64` path already established | None |
| `anndata` | No change — `ad.concat` + `ad.read_h5ad` already available | None |
| Plotly CDN tag | Fix stale version pin in `report_utils.py` line 517 | Code change only |

### What NOT to Add

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| Bokeh | Second interactive JS framework in same reports | Plotly (already present) |
| HoloViews / hvplot | Heavy deps, no gain over sc.pl.* + plotly for this use case | scanpy + plotly |
| Dash | Full app server — overkill for static HTML report generation | `fig.to_html(full_html=False)` |
| ipywidgets | Requires Jupyter kernel | base64 PNG + Plotly CDN |
| napari | No web-embeddable output | `sc.pl.spatial` → base64 PNG |
| vitessce | React-based viewer, requires build toolchain | Out of scope for CLI reports |
| `include_plotlyjs=True` inline | Adds ~3MB to every report; HPC reports are large already | CDN approach (existing) |
| `sq.pl.spatial_scatter` | squidpy's function is a thin wrapper over same matplotlib path | `sc.pl.spatial` (already used in codebase) |

### Python 3.10 Compatibility

`pyproject.toml` declares `requires-python = ">=3.10"`. All libraries relevant to v2.0:

| Package | Min Python | Notes |
|---------|-----------|-------|
| scanpy >=1.9 | 3.8+ | Confirmed via PyPI classifiers |
| matplotlib >=3.7 | 3.8+ | Confirmed |
| plotly >=5.18 | 3.8+ | Confirmed |
| anndata >=0.10 | 3.9+ | Confirmed |
| squidpy >=1.3 | 3.9+ | Confirmed |

No 3.11+ features are required for v2.0 plot generation or concat CLI. The existing codebase note
in the v1.0 section ("Python 3.11") reflected the sc_tools library consumers; the pyproject.toml
floor is 3.10 and all new v2.0 code must respect that floor.

### v2.0 Sources

- Codebase: `sc_tools/qc/report.py` lines 314–333 — existing `sc.pl.spatial` + `fig_to_base64` pattern (HIGH)
- Codebase: `sc_tools/qc/report_utils.py` lines 120–133, 517 — `plotly_to_html`, CDN injection (HIGH)
- Codebase: `sc_tools/pl/benchmarking.py` line 529 — existing `sc.pl.umap` pattern (HIGH)
- Codebase: `sc_tools/ingest/loaders.py` lines 903–947 — `concat_samples` implementation (HIGH)
- Codebase: `pyproject.toml` — `requires-python = ">=3.10"`, dep versions and extras (HIGH)
- [scanpy.pl.spatial docs](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.spatial.html) — `ax=`, `show=False` params stable since 1.8 (HIGH)
- [plotly.io.to_html docs](https://plotly.com/python-api-reference/generated/plotly.io.to_html.html) — `include_plotlyjs` options including `False` (CDN separate) (HIGH)

---
*v2.0 addendum: 2026-03-27*
