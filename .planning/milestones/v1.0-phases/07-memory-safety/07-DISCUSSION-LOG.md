# Phase 7: Memory Safety - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-03-24
**Phase:** 07-memory-safety
**Areas discussed:** IO Gateway tiers

---

## Area Selection

| Option | Description | Selected |
|--------|-------------|----------|
| IO Gateway tiers (Recommended) | How should the tiered loading work? Which operations get h5py-only vs backed vs full load? How does existing code factor in? | ✓ |
| Estimation model | How should `sct estimate` predict memory/runtime? Formula-based, empirical, or profile-based? | |
| Dry-run behavior | What should --dry-run report? Global flag or per-command? | |
| Scope & boundaries | Where does IO Gateway live? Replace all reads or CLI-only? | |

**User's choice:** IO Gateway tiers only
**Notes:** User was satisfied after 4 questions and chose to proceed to context creation without discussing remaining areas.

---

## IO Gateway Tiers

### Q1: Gateway location

| Option | Description | Selected |
|--------|-------------|----------|
| New sc_tools/io/ package (Recommended) | Dedicated module, clean separation from storage.py | ✓ |
| Extend storage.py | Add tiered loading to existing smart_read_h5ad() | |
| Extend memory/ package | Add gateway logic to existing memory/profiling.py | |

**User's choice:** New sc_tools/io/ package
**Notes:** None

### Q2: Gateway scope

| Option | Description | Selected |
|--------|-------------|----------|
| CLI-path only (Recommended) | Gateway used by sct CLI commands only. Library internals unchanged. | ✓ |
| All loading | Replace every anndata.read_h5ad call with gateway. | |
| CLI + MCP | Gateway for both CLI commands and MCP tool calls. | |

**User's choice:** CLI-path only
**Notes:** None

### Q3: Tier dispatch strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Operation-based (Recommended) | Each CLI command declares what access level it needs. Gateway reads declaration and loads accordingly. | ✓ |
| Size-based auto | Gateway checks file size on disk. Below threshold → full load, above → backed mode. | |
| Hybrid (size + operation) | Command declares minimum tier. Gateway promotes based on file size vs available memory. | |

**User's choice:** Operation-based
**Notes:** None

### Q4: Memory guard behavior

| Option | Description | Selected |
|--------|-------------|----------|
| Guard with override (Recommended) | Estimate before full load, refuse if exceeds available RAM, --force to override. | ✓ |
| Load as requested | Trust command's tier declaration, no pre-load check. | |
| Guard without override | Hard refuse, no --force escape hatch. | |

**User's choice:** Guard with override
**Notes:** None

---

## Claude's Discretion

- Estimation model design (formula-based vs empirical vs profile-based)
- --dry-run implementation strategy (global cli_handler interception vs per-command)
- --dry-run report contents (validation, planned operations, estimated resources)
- Backed mode (T2) interaction with partial X access
- Tier declaration mechanism (decorator parameter, enum, command metadata)
- Memory estimation formula calibration approach

## Deferred Ideas

None — discussion stayed within phase scope
