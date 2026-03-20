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
