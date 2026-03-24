# Phase 4: CLI Discovery - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md -- this log preserves the alternatives considered.

**Date:** 2026-03-23
**Phase:** 04-cli-discovery
**Areas discussed:** Schema format, Parameter metadata depth, Command naming convention, Output schema inclusion

---

## Schema format

| Option | Description | Selected |
|--------|-------------|----------|
| JSON Schema style | Each command's params as JSON Schema object. Standard format agents understand. | ✓ |
| OpenAPI-inspired | Model after OpenAPI operation objects. More verbose but familiar to web-tooling agents. | |
| Flat custom format | Simple custom schema with {name: {type, default, required, help}}. Minimal but non-standard. | |

**User's choice:** JSON Schema style
**Notes:** None -- straightforward selection of the recommended approach.

---

## Parameter metadata depth

| Option | Description | Selected |
|--------|-------------|----------|
| Standard | Type, default, required, description, enum. Auto-generated from Typer/Click. | ✓ |
| Rich with examples | Standard plus example values, semantic tags, constraints. Requires manual annotation. | |
| Minimal | Just name, type, required. Lean but agents need --help for the rest. | |

**User's choice:** Standard
**Notes:** Fully auto-generated from existing Typer parameter definitions -- no manual annotation burden.

---

## Command naming convention

| Option | Description | Selected |
|--------|-------------|----------|
| Space-separated | Match actual invocation: 'qc run', 'benchmark integration'. | ✓ |
| Dotted paths | Namespace-style: 'qc.run'. Requires translation for invocation. | |
| Flat slugs | Hyphenated: 'qc-run'. Single token but loses hierarchy. | |

**User's choice:** Space-separated
**Notes:** Consistent with how the CLI is actually invoked.

---

## Output schema inclusion

| Option | Description | Selected |
|--------|-------------|----------|
| Full contract | CLIResult output schema in $defs. Single source of truth for input + output. | ✓ |
| Input only | Only command parameters. Output is always CLIResult -- agent reads docs. | |
| You decide | Claude's discretion based on what Pydantic gives for free. | |

**User's choice:** Full contract
**Notes:** CLIResult schema auto-generated via model_json_schema(). sct schema becomes the single source of truth.

---

## Claude's Discretion

- Typer command tree walking strategy
- list-commands flat vs hierarchical
- describe command resolution (exact vs fuzzy)
- Top-level commands vs discovery subgroup
- Module placement (cli/discovery.py vs inline)

## Deferred Ideas

None -- discussion stayed within phase scope.
