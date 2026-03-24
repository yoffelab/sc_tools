# Phase 4: CLI Discovery - Research

**Researched:** 2026-03-23
**Domain:** Typer/Click introspection, JSON Schema generation, CLI metadata extraction
**Confidence:** HIGH

## Summary

Phase 4 adds three discovery commands (`list-commands`, `describe`, `schema`) that introspect the existing Typer command tree and Pydantic models to produce machine-readable catalogs. This is a metadata-only phase -- no pipeline logic, no heavy dependencies, no data loading.

The core technique is `typer.main.get_command(app)` which returns the underlying Click `Group` object. From there, `group.commands` gives all registered commands, `group.list_commands(ctx)` provides ordering, and each `Command.params` list gives full parameter metadata (name, type, default, required, help text, opts). Pydantic's `CLIResult.model_json_schema()` already produces the output schema with `$defs` for `Status`, `ErrorInfo`, and `Provenance`.

**Primary recommendation:** Create a single `sc_tools/cli/discovery.py` module with three commands registered as top-level commands on the main `app`. Walk the Click command tree recursively, map Click param types to JSON Schema types, and merge with `CLIResult.model_json_schema()` for the full contract. All three commands use `@cli_handler` for consistent CLIResult output.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- D-01: JSON Schema style. Each command's params described as a JSON Schema object (`type`, `default`, `description`, `enum`). Standard format agents already understand.
- D-02: `sct schema` emits a top-level object with `commands` keyed by space-separated name, each containing a `params` JSON Schema + output schema. Shared types go in `$defs`.
- D-03: Standard depth: type, default, required, description, and enum values where applicable. Fully auto-generated from Typer/Click internals -- no manual annotation needed.
- D-04: No rich extras (examples, semantic tags, constraints beyond enum). Keep it maintainable by deriving everything from the existing Typer parameter definitions.
- D-05: Space-separated names matching actual CLI invocation: `"qc run"`, `"benchmark integration"`, `"validate"`, `"status"`. Agent can split on space to get group + subcommand.
- D-06: Top-level commands (no subcommand) use bare name: `"validate"`, `"status"`, `"report"`.
- D-07: Full contract -- include CLIResult output schema as a shared `$defs` entry. `sct schema` is the single source of truth for both input params and output shape.
- D-08: CLIResult schema auto-generated via `CLIResult.model_json_schema()` from Pydantic. Per-command `data` field shape documented where it varies.

### Claude's Discretion
- How to walk the Typer command tree (Click's internal `commands` dict, `get_params()`, etc.)
- Whether `list-commands` is a flat list or includes group hierarchy
- How `describe` resolves command names (exact match vs fuzzy)
- Whether to add these as top-level commands or a `discovery` subgroup
- Internal module placement (new `cli/discovery.py` or inline in `cli/__init__.py`)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| DSC-01 | `sct list-commands --json` -- machine-readable catalog of all commands with params, types, defaults | Click tree walk via `typer.main.get_command(app)` provides all metadata; verified on actual app |
| DSC-02 | `sct describe <cmd>` -- JSON schema for specific command's params and output format | Click `Command.params` + `CLIResult.model_json_schema()` combine to produce per-command schema |
| DSC-03 | `sct schema` -- full CLI contract as JSON document (Typer command tree + Pydantic output schemas) | Merge of walk-all-commands + Pydantic schema with `$defs`; verified both APIs produce correct output |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- Output paths: active projects live at `~/Documents/projects/active/<project>/` -- never repo root
- Phase transitions: call `set_phase_status` on completion, dispatch `documentor` to update Plan.md and Journal.md
- TDD skill (Tier 1): No production code without a failing test first; red-green-refactor cycle mandatory
- Linting: Ruff; never commit failing lint
- Lazy imports (CLI-08): Heavy dependencies loaded at command execution, not startup
- Test files: `sc_tools/tests/test_{module}.py`

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| typer | 0.24.1 | CLI framework (already installed) | Project standard since Phase 2 |
| click | 8.3.0 | Underlying CLI library (Typer wraps it) | Used for introspection APIs |
| pydantic | 2.x | CLIResult model + JSON Schema generation | Project standard since Phase 2 |

### Supporting
No new dependencies required. Everything is available from the existing stack.

**Installation:** None required -- all dependencies already installed.

## Architecture Patterns

### Recommended Project Structure
```
sc_tools/cli/
    __init__.py          # existing -- add import + registration for discovery
    discovery.py         # NEW -- list-commands, describe, schema commands
    qc.py               # existing
    preprocess.py        # existing
    validate.py          # existing
    benchmark.py         # existing
    status.py            # existing
sc_tools/tests/
    test_cli_discovery.py  # NEW -- tests for discovery commands
```

### Pattern 1: Click Tree Walk
**What:** Convert Typer app to Click Group, then walk recursively to extract command metadata.
**When to use:** All three discovery commands need this.
**Example:**
```python
# Source: verified on actual sc_tools app (see research introspection output)
import click
import typer

def _walk_commands(app: typer.Typer) -> dict[str, click.Command]:
    """Walk the Typer command tree, return flat dict of space-separated name -> Click Command."""
    click_app = typer.main.get_command(app)
    commands: dict[str, click.Command] = {}
    ctx = click.Context(click_app)

    for name in click_app.list_commands(ctx):
        cmd = click_app.commands[name]
        if isinstance(cmd, click.Group):
            grp_ctx = click.Context(cmd, parent=ctx)
            for sub_name in cmd.list_commands(grp_ctx):
                sub_cmd = cmd.get_command(grp_ctx, sub_name)
                commands[f"{name} {sub_name}"] = sub_cmd
        else:
            commands[name] = cmd

    return commands
```

### Pattern 2: Click Param to JSON Schema
**What:** Map Click parameter objects to JSON Schema property definitions.
**When to use:** Building the `params` schema for each command.
**Example:**
```python
# Source: verified against actual Click param objects from sc_tools app
_CLICK_TYPE_MAP = {
    "text": "string",
    "integer": "integer",
    "float": "number",
    "boolean": "boolean",
    "path": "string",
    "uuid": "string",
    "choice": "string",
}

def _param_to_schema(param: click.Parameter) -> dict:
    """Convert a Click parameter to a JSON Schema property."""
    schema: dict = {}
    type_name = param.type.name
    schema["type"] = _CLICK_TYPE_MAP.get(type_name, "string")

    if param.help:
        schema["description"] = param.help
    if param.default is not None:
        schema["default"] = param.default
    if hasattr(param.type, "choices") and param.type.choices:
        schema["enum"] = list(param.type.choices)

    return schema
```

### Pattern 3: Consistent Registration as Top-Level Commands
**What:** Register discovery commands directly on the main `app` (not a subgroup).
**When to use:** These are utility commands like `version`, not pipeline commands.
**Example:**
```python
# In cli/discovery.py:
# Register on the main app, following the same pattern as version command

# In cli/__init__.py (at bottom with other imports):
from sc_tools.cli.discovery import discovery_app  # noqa: E402
# Or register individual commands directly
```

### Discretion Decisions (Recommendations)

1. **Tree walk method:** Use `typer.main.get_command(app)` to get Click Group, then `commands` dict + `list_commands()`. Verified working on actual app. (HIGH confidence)

2. **list-commands format:** Flat list with group info as metadata. Each entry has `name` (space-separated), `group` (or null for top-level), and `params`. This is simpler for agents to consume than nested hierarchy. (MEDIUM confidence -- either works, flat is simpler)

3. **describe resolution:** Exact match only. Space-separated name (`"qc run"`) must match exactly. On mismatch, return an error CLIResult with available command names as suggestion. No fuzzy matching -- agents should use `list-commands` first. (HIGH confidence)

4. **Registration:** Top-level commands on `app` (not a `discovery` subgroup). Three commands is not enough to warrant a subgroup, and `sct list-commands` reads better than `sct discovery list-commands`. (MEDIUM confidence)

5. **Module placement:** New `cli/discovery.py` file (not inline in `__init__.py`). Follows established pattern where each command domain has its own module. (HIGH confidence)

6. **Exclude discovery commands from their own output:** The three discovery commands (`list-commands`, `describe`, `schema`) should NOT appear in the catalog. They are meta-commands, not pipeline operations. Agents discovering the CLI want to know what operations they can perform, not how to discover operations. (MEDIUM confidence -- could go either way, but excluding is cleaner)

### Anti-Patterns to Avoid
- **Parsing help text:** Never parse `--help` output. Use Click's structured param objects directly.
- **Manual command registry:** Don't maintain a separate list of commands. Walk the Typer tree at runtime so it's always in sync.
- **Heavy imports in discovery:** Discovery commands must not import scanpy, torch, etc. They only introspect the CLI structure.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| JSON Schema for CLIResult | Manual schema dict | `CLIResult.model_json_schema()` | Pydantic v2 auto-generates correct JSON Schema with `$defs` |
| Command tree enumeration | Manual command list | `typer.main.get_command(app).commands` | Always in sync with actual registered commands |
| Parameter type mapping | Complex type introspection | `param.type.name` + simple map | Click normalizes all types to a small set of named types |
| Choice/enum extraction | Parsing help text | `param.type.choices` | Click Choice type exposes choices directly |

**Key insight:** Both Typer/Click and Pydantic already expose structured metadata. The discovery commands are thin wrappers that format existing metadata as JSON Schema.

## Common Pitfalls

### Pitfall 1: Circular Import
**What goes wrong:** `discovery.py` imports from `cli/__init__.py` which imports from `discovery.py`.
**Why it happens:** Discovery needs the `app` object to introspect, and `__init__.py` needs to register discovery commands.
**How to avoid:** Discovery module defines its commands. `__init__.py` imports and registers them at the bottom (same pattern as qc, benchmark, etc.). Discovery functions receive `app` as parameter or import it lazily inside the command function body.
**Warning signs:** ImportError at module load time.

### Pitfall 2: Group Callbacks Appearing as Commands
**What goes wrong:** `qc_callback`, `status_callback` etc. show up as commands when walking the tree.
**Why it happens:** Typer registers callbacks differently than commands, but Click tree may still expose them.
**How to avoid:** Filter commands -- skip entries where the command is a Group (check `isinstance(cmd, click.Group)`). Only include leaf commands that have actual functionality.
**Warning signs:** Commands like `"qc"` (with no subcommand) appearing alongside `"qc run"`.

### Pitfall 3: `--help` and `--version` Params in Output
**What goes wrong:** Click auto-adds `--help` param to every command. `--version` is on the main app callback.
**Why it happens:** Click's built-in params.
**How to avoid:** Filter out params named `help` and `version` (or any param with `is_eager=True` and `expose_value=False`). These are Click internals, not command parameters.
**Warning signs:** Every command showing `help: {type: boolean, default: false}` in its schema.

### Pitfall 4: validate_run Bypasses cli_handler
**What goes wrong:** `validate run` doesn't use `@cli_handler` -- it manages `_emit` and `SystemExit` directly. Its behavior is the same from the outside but its function doesn't return CLIResult.
**Why it happens:** Phase 3 decision -- custom exit code 2 for data errors.
**How to avoid:** Discovery commands introspect Click params (which are identical regardless of handler pattern). No special handling needed -- the introspection layer doesn't care about the handler decorator.
**Warning signs:** None for discovery; this is informational context.

### Pitfall 5: Stub Groups (integrate, celltype) Have No Subcommands
**What goes wrong:** Walking `integrate` and `celltype` groups finds zero commands (they're stubs from Phase 2).
**Why it happens:** No Phase 3 commands were created for these groups.
**How to avoid:** Include the group in the catalog with an empty `commands` list, or skip groups with no leaf commands. Recommend: skip them entirely -- agents should only see actionable commands.
**Warning signs:** Empty groups cluttering the catalog.

## Code Examples

Verified patterns from actual codebase introspection:

### Walking the Command Tree (verified)
```python
# Tested against sc_tools app -- produces correct output
import click
import typer
from sc_tools.cli import app

click_app = typer.main.get_command(app)
ctx = click.Context(click_app)

# Result: ['version', 'qc', 'report', 'preprocess', 'validate',
#          'status', 'integrate', 'benchmark', 'celltype']
names = click_app.list_commands(ctx)

# Subcommands within groups:
# qc -> run
# report -> generate
# preprocess -> run
# validate -> run
# status -> show
# benchmark -> integration
# integrate -> (empty)
# celltype -> (empty)
# version -> (top-level, not a group)
```

### Extracting Parameter Metadata (verified)
```python
# Actual output from qc run command params:
# file: Argument, type=text, required=True, help="Path to h5ad checkpoint file"
# project_dir: Option, type=text, required=False, default=".", opts=["--project-dir"]
# modality: Option, type=text, required=False, default="visium", opts=["--modality", "-m"]
# sample_col: Option, type=text, required=False, default="library_id", opts=["--sample-col"]
```

### CLIResult JSON Schema (verified)
```python
# CLIResult.model_json_schema() produces:
# {
#   "$defs": {
#     "ErrorInfo": { ... },
#     "Provenance": { ... },
#     "Status": { "enum": ["success", "error", "partial", "skipped"], "type": "string" }
#   },
#   "properties": {
#     "status": { "$ref": "#/$defs/Status" },
#     "command": { "type": "string" },
#     "data": { "type": "object" },
#     "artifacts": { "items": { "type": "string" }, "type": "array" },
#     "provenance": { "$ref": "#/$defs/Provenance" },
#     "message": { "type": "string" },
#     "error": { ... },
#     "partial_failures": { ... }
#   },
#   "required": ["status", "command", "provenance"]
# }
```

### Complete list-commands Expected Output Shape
```json
{
  "status": "success",
  "command": "list-commands",
  "data": {
    "commands": [
      {
        "name": "qc run",
        "group": "qc",
        "description": "Run QC metrics on a checkpoint file, output JSON summary",
        "params": {
          "type": "object",
          "properties": {
            "file": { "type": "string", "description": "Path to h5ad checkpoint file" },
            "project_dir": { "type": "string", "description": "Project directory", "default": "." },
            "modality": { "type": "string", "description": "Data modality ...", "default": "visium" },
            "sample_col": { "type": "string", "description": "Sample column in obs", "default": "library_id" }
          },
          "required": ["file"]
        }
      }
    ],
    "count": 7
  },
  "provenance": { "command": "list-commands", ... },
  "message": "7 commands available"
}
```

### Complete schema Expected Output Shape
```json
{
  "status": "success",
  "command": "schema",
  "data": {
    "schema_version": "1.0",
    "sc_tools_version": "...",
    "commands": {
      "qc run": {
        "group": "qc",
        "description": "...",
        "params": { "type": "object", "properties": { ... }, "required": [ ... ] }
      },
      "benchmark integration": { ... }
    },
    "output": {
      "$ref": "#/$defs/CLIResult"
    },
    "$defs": {
      "CLIResult": { ... },
      "Status": { ... },
      "ErrorInfo": { ... },
      "Provenance": { ... }
    }
  },
  "provenance": { ... },
  "message": "Full CLI schema: 7 commands, 4 shared types"
}
```

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest (already configured) |
| Config file | `pyproject.toml` (pytest section) |
| Quick run command | `pytest sc_tools/tests/test_cli_discovery.py -v --tb=short -q` |
| Full suite command | `pytest sc_tools/tests/ -v --tb=short -q` |

### Phase Requirements to Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| DSC-01 | `list-commands --json` returns catalog with all commands, params, types, defaults | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestListCommands -x` | Wave 0 |
| DSC-01 | Catalog includes correct param types and defaults for known commands | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestListCommands::test_param_types -x` | Wave 0 |
| DSC-02 | `describe qc run` returns JSON schema for that command | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestDescribe -x` | Wave 0 |
| DSC-02 | `describe nonexistent` returns error CLIResult with suggestion | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestDescribe::test_unknown_command -x` | Wave 0 |
| DSC-03 | `schema` returns full contract with commands + $defs | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestSchema -x` | Wave 0 |
| DSC-03 | Schema $defs includes CLIResult, Status, ErrorInfo, Provenance | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestSchema::test_defs -x` | Wave 0 |
| ALL | Discovery commands follow cli_handler pattern (CLIResult envelope) | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestOutputFormat -x` | Wave 0 |
| ALL | No heavy imports (scanpy, torch) in discovery module | unit | `pytest sc_tools/tests/test_cli_discovery.py::TestLazyImports -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_cli_discovery.py -v --tb=short -q`
- **Per wave merge:** `pytest sc_tools/tests/ -v --tb=short -q`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/test_cli_discovery.py` -- covers DSC-01, DSC-02, DSC-03
- No framework install needed -- pytest already configured
- No new fixtures needed -- uses `CliRunner` from `typer.testing` (already in test_cli.py)

## Sources

### Primary (HIGH confidence)
- Actual codebase introspection -- `python3` scripts run against `sc_tools.cli.app` to verify Click tree walk, param extraction, and type mapping
- `sc_tools/cli/__init__.py` -- Typer app structure, command registration pattern, cli_handler decorator
- `sc_tools/models/result.py` -- CLIResult, Status, ErrorInfo, Provenance models + `model_json_schema()` output
- `sc_tools/cli/qc.py`, `benchmark.py`, `preprocess.py`, `validate.py`, `status.py` -- all command implementations (introspection targets)
- `sc_tools/tests/test_cli.py` -- existing test patterns for CLI testing with CliRunner

### Secondary (MEDIUM confidence)
- [Typer commands tutorial](https://typer.tiangolo.com/tutorial/commands/) -- command registration patterns
- [Click Groups documentation](https://click.palletsprojects.com/en/stable/commands/) -- Group.list_commands, Command.params API

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - no new dependencies, everything verified installed and working
- Architecture: HIGH - all introspection APIs verified against actual codebase
- Pitfalls: HIGH - tested actual tree walk, identified real edge cases (stub groups, callbacks, help params)

**Research date:** 2026-03-23
**Valid until:** 2026-04-23 (stable -- Typer/Click APIs unlikely to change)
