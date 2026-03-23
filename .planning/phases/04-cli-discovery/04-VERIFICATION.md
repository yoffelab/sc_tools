---
phase: 04-cli-discovery
verified: 2026-03-23T00:00:00Z
status: passed
score: 5/5 must-haves verified
re_verification: false
---

# Phase 4: CLI Discovery Verification Report

**Phase Goal:** Agents can programmatically discover all available commands, their parameters, and output schemas without parsing help text
**Verified:** 2026-03-23
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `sct list-commands --json` returns a JSON catalog with every leaf command's name, params, types, defaults | VERIFIED | Returns 6 commands (benchmark integration, preprocess run, qc run, report generate, status show, validate run) with full JSON Schema params. Confirmed by TestListCommands (6 tests, all pass) and spot-check: count=6, names correct. |
| 2 | `sct describe <cmd>` returns JSON schema for a specific command's input params and output format | VERIFIED | `sct describe "qc run"` returns status=success, data.name="qc run", data.params with properties, data.output_schema with $defs (Status, ErrorInfo, Provenance). TestDescribe (3 tests, all pass). |
| 3 | `sct schema` returns the full CLI contract as a single JSON document with command tree and Pydantic $defs | VERIFIED | Returns commands dict keyed by space-separated names, $defs with CLIResult, Status, ErrorInfo, Provenance, schema_version="1.0", sc_tools_version non-empty. TestSchema (5 tests, all pass). |
| 4 | Discovery commands produce valid CLIResult JSON envelope (same as all other commands) | VERIFIED | All three commands return JSON with status, command, provenance keys. TestOutputFormat passes. Spot-check confirmed all three commands exit 0 with valid CLIResult. |
| 5 | Discovery module imports no heavy dependencies (scanpy, torch, scvi-tools) | VERIFIED | No `import scanpy`, `import torch`, `import scvi` found in discovery.py. TestLazyImports passes — checks sys.modules before/after import. |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `sc_tools/cli/discovery.py` | list-commands, describe, schema command implementations + _walk_commands + _param_to_schema helpers; min 80 lines | VERIFIED | 232 lines. Contains _walk_commands, _param_to_schema, _command_to_entry, register_discovery with all three commands. No heavy imports. |
| `sc_tools/tests/test_cli_discovery.py` | Tests for all three discovery commands; min 80 lines | VERIFIED | 170 lines. Contains TestListCommands (6 tests), TestDescribe (3 tests), TestSchema (5 tests), TestOutputFormat (1 test), TestLazyImports (1 test) = 16 tests total, all pass. |
| `sc_tools/cli/__init__.py` | Registration of discovery commands on main app; contains `from sc_tools.cli.discovery` | VERIFIED | Line 274: `from sc_tools.cli.discovery import register_discovery` followed by `register_discovery(app)` on line 275. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `sc_tools/cli/discovery.py` | `sc_tools/cli/__init__.py` | `register_discovery(app)` called at bottom of __init__.py after app is fully constructed | WIRED | Line 274-275 of __init__.py: `from sc_tools.cli.discovery import register_discovery` then `register_discovery(app)`. Circular import avoided — app is imported lazily inside each command body via `from sc_tools.cli import app as main_app`. |
| `sc_tools/cli/discovery.py` | `sc_tools/models/result.py` | Uses CLIResult, Status, Provenance for output + model_json_schema() for $defs | WIRED | discovery.py line 18: `from sc_tools.models.result import CLIResult, Provenance, Status, _get_version`. `CLIResult.model_json_schema()` called on lines 187 and 211. |
| `sc_tools/tests/test_cli_discovery.py` | `sc_tools/cli/__init__.py` | CliRunner invokes app commands | WIRED | Line 15: `from sc_tools.cli import app`. Module-level `runner = CliRunner()`. All tests use `runner.invoke(app, [...])`. |

### Data-Flow Trace (Level 4)

Discovery commands are introspection tools, not data-rendering components — they read the live Typer/Click command tree at invocation time rather than a database or external store. Level 4 trace is still meaningful to confirm the introspection is connected to real command metadata rather than hardcoded output.

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|--------------|--------|--------------------|--------|
| `discovery.py` — `list_commands` | `commands` dict from `_walk_commands(main_app)` | `typer.main.get_command(target_app)` — walks live Click tree at runtime | Yes — 6 real commands discovered from registered Typer app; `qc run` params include file=string, modality=visium default, all from actual Click param definitions | FLOWING |
| `discovery.py` — `describe` | `entry` dict from `_command_to_entry` + `CLIResult.model_json_schema()` | Same live Click tree + Pydantic introspection | Yes — output_schema $defs contain Status, ErrorInfo, Provenance pulled from actual Pydantic model | FLOWING |
| `discovery.py` — `schema` | `commands_dict` + `defs` from `CLIResult.model_json_schema()` | Live Click tree + Pydantic model | Yes — $defs keys confirmed as CLIResult, ErrorInfo, Provenance, Status by spot-check | FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| `sct list-commands --json` returns 6 real commands | `python -c "... r.invoke(app, ['list-commands', '--json'])"` | count=6, names=['benchmark integration', 'preprocess run', 'qc run', 'report generate', 'status show', 'validate run'] | PASS |
| `sct describe "qc run"` output_schema has $defs | `python -c "... '$defs' in d2['data']['output_schema']"` | True | PASS |
| `sct schema` $defs keys are correct | `python -c "... sorted(d3['data']['$defs'].keys())"` | ['CLIResult', 'ErrorInfo', 'Provenance', 'Status'] | PASS |
| All 16 discovery tests pass | `python -m pytest sc_tools/tests/test_cli_discovery.py -v` | 16 passed in 0.33s | PASS |
| All 24 existing CLI tests still pass | `python -m pytest sc_tools/tests/test_cli.py -v` | 24 passed in 0.35s | PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| DSC-01 | 04-01-PLAN.md | `sct list-commands --json` — machine-readable catalog of all commands with params, types, defaults | SATISFIED | `list-commands` command in discovery.py, returns 6 commands with full JSON Schema params including types and defaults. TestListCommands 6/6 pass. |
| DSC-02 | 04-01-PLAN.md | `sct describe <cmd>` — JSON schema for specific command's params and output format | SATISFIED | `describe` command in discovery.py, returns per-command JSON schema plus CLIResult output_schema with $defs. TestDescribe 3/3 pass. |
| DSC-03 | 04-01-PLAN.md | `sct schema` — full CLI contract as JSON document (Typer command tree + Pydantic output schemas) | SATISFIED | `schema` command in discovery.py, returns full contract with commands dict keyed by space-separated name, $defs with CLIResult/Status/ErrorInfo/Provenance, schema_version, sc_tools_version. TestSchema 5/5 pass. |

**No orphaned requirements for Phase 4.** REQUIREMENTS.md traceability table maps exactly DSC-01, DSC-02, DSC-03 to Phase 4 — same as PLAN frontmatter. All three marked Complete.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| — | — | None found | — | — |

No TODOs, FIXMEs, placeholders, heavy imports, empty returns, or hardcoded stubs detected in discovery.py or test_cli_discovery.py.

### Human Verification Required

None. All behaviors are programmatically verifiable via CliRunner invocations and module introspection. No UI, real-time, or external service components exist in this phase.

### Gaps Summary

No gaps. All five must-have truths are fully verified:

1. All three discovery commands exist and are substantive (232-line module, no stubs).
2. All commands are wired — registered on the live app via `register_discovery(app)` at the bottom of `__init__.py`.
3. Data flows from the real Typer/Click command tree — 6 commands discovered, param metadata is accurate (confirmed by asserting `modality` default is `"visium"` from the actual `qc run` Click definition).
4. CLIResult envelope is used by all three commands, with correct exit codes (0 success, 1 user error).
5. No heavy imports — lazy import pattern used correctly, `app` imported inside function bodies to avoid circular dependency.

Phase goal is fully achieved: agents can call `sct list-commands --json`, `sct describe <cmd>`, and `sct schema` to obtain complete machine-readable CLI metadata without parsing help text.

---

_Verified: 2026-03-23_
_Verifier: Claude (gsd-verifier)_
