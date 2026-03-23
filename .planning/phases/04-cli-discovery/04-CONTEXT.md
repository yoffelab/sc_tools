# Phase 4: CLI Discovery - Context

**Gathered:** 2026-03-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Agents can programmatically discover all available commands, their parameters, and output schemas without parsing help text. Requirements DSC-01, DSC-02, DSC-03. Three new commands: `sct list-commands`, `sct describe`, `sct schema`. These introspect the existing Typer command tree and Pydantic models from Phases 2-3 — no new pipeline logic.

</domain>

<decisions>
## Implementation Decisions

### Schema format
- **D-01:** JSON Schema style. Each command's params described as a JSON Schema object (`type`, `default`, `description`, `enum`). Standard format agents already understand.
- **D-02:** `sct schema` emits a top-level object with `commands` keyed by space-separated name, each containing a `params` JSON Schema + output schema. Shared types go in `$defs`.

### Parameter metadata depth
- **D-03:** Standard depth: type, default, required, description, and enum values where applicable. Fully auto-generated from Typer/Click internals — no manual annotation needed.
- **D-04:** No rich extras (examples, semantic tags, constraints beyond enum). Keep it maintainable by deriving everything from the existing Typer parameter definitions.

### Command naming convention
- **D-05:** Space-separated names matching actual CLI invocation: `"qc run"`, `"benchmark integration"`, `"validate"`, `"status"`. Agent can split on space to get group + subcommand.
- **D-06:** Top-level commands (no subcommand) use bare name: `"validate"`, `"status"`, `"report"`.

### Output schema inclusion
- **D-07:** Full contract — include CLIResult output schema as a shared `$defs` entry. `sct schema` is the single source of truth for both input params and output shape.
- **D-08:** CLIResult schema auto-generated via `CLIResult.model_json_schema()` from Pydantic. Per-command `data` field shape documented where it varies (e.g., qc metrics dict vs benchmark rankings).

### Claude's Discretion
- How to walk the Typer command tree (Click's internal `commands` dict, `get_params()`, etc.)
- Whether `list-commands` is a flat list or includes group hierarchy
- How `describe` resolves command names (exact match vs fuzzy)
- Whether to add these as top-level commands or a `discovery` subgroup
- Internal module placement (new `cli/discovery.py` or inline in `cli/__init__.py`)

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### CLI foundation (Phase 2-3 output)
- `sc_tools/cli/__init__.py` -- Typer app, command group registration, `cli_handler`, `_emit`, `_render_rich`
- `sc_tools/models/result.py` -- CLIResult, Status, ErrorInfo, Provenance Pydantic models (use `model_json_schema()`)
- `sc_tools/errors.py` -- Exception hierarchy

### Command implementations (introspection targets)
- `sc_tools/cli/qc.py` -- `sct qc run`, `sct report` commands
- `sc_tools/cli/preprocess.py` -- `sct preprocess run`
- `sc_tools/cli/benchmark.py` -- `sct benchmark integration`
- `sc_tools/cli/validate.py` -- `sct validate`
- `sc_tools/cli/status.py` -- `sct status`

### Pipeline metadata
- `sc_tools/pipeline.py` -- `STANDARD_PHASES`, `PhaseSpec` objects with checkpoint paths and dependencies

### Requirements
- `.planning/REQUIREMENTS.md` -- DSC-01, DSC-02, DSC-03

### Prior context
- `.planning/phases/02-cli-foundation/02-CONTEXT.md` -- CLIResult design, error taxonomy
- `.planning/phases/03-core-commands/03-CONTEXT.md` -- Command structure, cli/ package layout

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `CLIResult.model_json_schema()` -- Pydantic v2 auto-generates JSON Schema for the output envelope
- Typer's Click internals: `app.registered_commands`, Click `Command.params`, `Command.help` -- all param metadata is already there
- `_get_version()` in `models/result.py` -- version string for schema metadata

### Established Patterns
- Lazy imports inside command functions (CLI-08) -- discovery commands must also follow this
- `cli_handler` decorator on every command function -- discovery commands should use it too for consistent CLIResult output
- `add_typer()` registration in `__init__.py` -- 8 groups already registered

### Integration Points
- New discovery commands register on `app` directly (top-level) or via a new Typer group
- Typer command tree is the single source of truth -- introspect `app` object at runtime
- Pydantic `model_json_schema()` provides `$defs` for CLIResult, ErrorInfo, Provenance, Status

</code_context>

<specifics>
## Specific Ideas

No specific requirements -- open to standard approaches for Typer introspection.

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 04-cli-discovery*
*Context gathered: 2026-03-23*
