# Phase 5: Provenance - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md -- this log preserves the alternatives considered.

**Date:** 2026-03-23
**Phase:** 05-provenance
**Areas discussed:** Sidecar writing strategy, Provenance model fields, Lineage trace design, Leiden reproducibility, File portability

---

## Sidecar Writing Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Auto via cli_handler | After successful exec, write .provenance.json next to each artifact. No opt-out. | ✓ |
| Explicit per-command | Each command writes its own sidecar with custom fields. More boilerplate. | |
| Auto with --no-provenance opt-out | Same as auto but with opt-out flag for exploratory runs. | |

**User's choice:** Auto via cli_handler
**Notes:** Minimal code change — cli_handler already has access to CLIResult and artifacts list.

---

## Provenance Model Fields (Input Files)

| Option | Description | Selected |
|--------|-------------|----------|
| Explicit params only | Record inputs from CLI flags. SHA256 at read time. | |
| Hybrid: sidecar + adata.uns | Always write sidecar. Also embed in adata.uns for h5ad outputs. | ✓ |
| Minimal (no checksums) | Record paths but skip SHA256. | |

**User's choice:** Hybrid (sidecar + adata.uns)
**Notes:** User raised concern about file portability across storage locations. Hybrid ensures provenance survives file moves — sidecar is canonical, adata.uns is portable backup. This was a user-initiated gray area not in the original list.

---

## Provenance Param Serialization

| Option | Description | Selected |
|--------|-------------|----------|
| All flags as-passed | Every CLI flag including defaults. Agent can replay exact command. | ✓ |
| Non-default only | Only explicitly set flags. Smaller but loses default context. | |
| All flags + environment | Flags plus env vars (CUDA, SCT_*). Most complete but noisy. | |

**User's choice:** All flags as-passed
**Notes:** None.

---

## Lineage Trace Design

| Option | Description | Selected |
|--------|-------------|----------|
| Follow sidecar input refs | Recurse through sidecar inputs. Flat chronological list output. | ✓ |
| Graph with DOT export | Full DAG with JSON + DOT visualization. | |
| You decide | Claude's discretion. | |

**User's choice:** Follow sidecar input refs
**Notes:** Flat list ordered chronologically. Stop at raw data origins (no sidecar).

---

## Leiden Reproducibility (PRV-05)

| Option | Description | Selected |
|--------|-------------|----------|
| Thread through call stack | Add random_state param to all Leiden-calling functions. Default 0. | ✓ |
| Global seed at CLI boundary | Set numpy/random seed once. Simpler but fragile. | |
| You decide | Claude's discretion. | |

**User's choice:** Thread through call stack
**Notes:** Record resolution and random_state in provenance. Test identical params = identical results.

---

## File Portability (User-Initiated)

| Option | Description | Selected |
|--------|-------------|----------|
| Relative paths + SHA256 matching | Relative paths in sidecars. SHA256 fallback when path missing. | ✓ |
| Absolute paths only | Simple but breaks on file moves. | |
| Content-addressed (SHA256 key) | Most robust but requires file indexing. Overkill. | |

**User's choice:** Relative paths + SHA256 matching
**Notes:** User has storage issues requiring file moves. Trace falls back to: (1) adata.uns, (2) SHA256 scan, (3) mark as origin.

---

## Claude's Discretion

- SHA256 computation strategy for large h5ad files
- Provenance model internal structure (nested vs flat)
- How cli_handler passes input file info to provenance
- adata.uns provenance structure
- Trace depth limit / cycle detection
- peak_memory_mb measurement approach

## Deferred Ideas

- W3C PROV-JSON export (v2 INF-01)
- RO-Crate packaging (v2 INF-02)
- DB-backed provenance (v2 INF-03)
- DOT/graph lineage visualization
