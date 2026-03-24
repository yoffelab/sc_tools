---
phase: 08-multi-omic-assembly
plan: 02
subsystem: assembly-embed-cli
tags: [embedding, mofa, multivi, totalvi, cli, mudata, multi-omic]
dependency_graph:
  requires: [08-01]
  provides: [embedding-backends, assemble-cli]
  affects: [sc_tools/assembly/embed/, sc_tools/cli/assemble.py]
tech_stack:
  added: [mofapy2, muon-embedding]
  patterns: [EmbeddingBackend-Protocol, register-pattern, cli-handler]
key_files:
  created:
    - sc_tools/assembly/embed/__init__.py
    - sc_tools/assembly/embed/_base.py
    - sc_tools/assembly/embed/_mofa.py
    - sc_tools/assembly/embed/_multivi.py
    - sc_tools/assembly/embed/_totalvi.py
    - sc_tools/cli/assemble.py
    - sc_tools/tests/test_assembly_embed.py
    - sc_tools/tests/test_cli_assemble.py
  modified:
    - sc_tools/assembly/_atlas.py
    - sc_tools/assembly/__init__.py
    - sc_tools/cli/__init__.py
decisions:
  - "Return type -> None on CLI commands to avoid Typer get_type_hints NameError on Python 3.10"
  - "MOFA+ integration tests guarded by both muon and mofapy2 importorskip"
  - "register_assemble(app) pattern with nested Typer subapp matches existing CLI registration"
metrics:
  duration: 6min
  completed: "2026-03-24T22:50:00Z"
  tasks: 2
  files: 11
  tests_added: 14
  tests_passing: 22
  tests_skipped: 3
---

# Phase 08 Plan 02: Embedding Backends + CLI Summary

Joint embedding backends (MOFA+, MultiVI, TotalVI) with pluggable dispatch registry and `sct assemble` CLI command group exposing build/embed/query subcommands.

## One-liner

EmbeddingBackend Protocol with MOFA+/MultiVI/TotalVI dispatch, plus `sct assemble build|embed|query` CLI wrapping MultiOmicAtlas operations.

## What Was Built

### Task 1: Embedding Backend Protocol, Registry, and Three Backends

- **EmbeddingBackend Protocol** (`_base.py`): `runtime_checkable` Protocol with `run(mdata, *, n_factors=15, **kwargs) -> tuple[ndarray, dict]`, mirroring `CelltypeBackend` pattern (D-08)
- **Registry**: `register_embedding_backend`, `get_embedding_backend` (raises ValueError with available methods), `list_embedding_methods` (sorted)
- **MofaBackend** (`_mofa.py`): Calls `muon.tl.mofa(mdata, use_obs="union", n_factors=n_factors)` (D-01 outer join, D-07 recommended default)
- **MultiviBackend** (`_multivi.py`): Validates RNA+ATAC modalities present, raises `SCToolsDataError` with clear message if not (Pitfall 3)
- **TotalviBackend** (`_totalvi.py`): Validates RNA+protein modalities present, same error pattern
- **MultiOmicAtlas.embed()**: New method dispatching to backend by name string
- **8 tests**: 6 pass, 2 skip (mofapy2 not installed)

### Task 2: sct assemble CLI Command Group

- **`sct assemble build`**: Takes h5ad files with `--modality` flags, builds h5mu via `MultiOmicAtlas.from_modalities`, validates input count matches modality count
- **`sct assemble embed`**: Loads h5mu, runs embedding via `atlas.embed(method=...)`, saves updated atlas
- **`sct assemble query`**: Loads h5mu, returns celltype proportions as JSON, optional `--patient` filter
- **Registration**: `register_assemble(app)` pattern in `cli/__init__.py`, `mudata`/`muon`/`mofapy2` added to `_DEP_INSTALL`
- **6 CLI tests**: 5 pass, 1 skip (mofapy2 not installed)

## Commits

| Task | Commit | Description |
|------|--------|-------------|
| 1 (RED) | 731c081 | Failing tests for embedding backends |
| 1 (GREEN) | 5404cfe | Embedding backend Protocol, registry, three backends |
| 2 | 472798e | sct assemble CLI command group (build, embed, query) |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Return type annotation causes Typer NameError on Python 3.10**
- **Found during:** Task 2
- **Issue:** `-> CLIResult` return type annotation on CLI functions causes `NameError: name 'CLIResult' is not defined` when Typer evaluates type hints at import time (known Phase 03 issue)
- **Fix:** Changed return type to `-> None` on all three CLI command functions
- **Files modified:** sc_tools/cli/assemble.py

**2. [Rule 1 - Bug] MOFA+ tests need mofapy2 guard in addition to muon**
- **Found during:** Task 1 (GREEN phase)
- **Issue:** `muon` is installed but `mofapy2` is not, causing ImportError at MOFA+ runtime
- **Fix:** Added `pytest.importorskip("mofapy2")` alongside muon guard for MOFA+ integration tests
- **Files modified:** sc_tools/tests/test_assembly_embed.py

## Known Stubs

None. All code paths are functional (MOFA+ integration deferred to environments with mofapy2 installed).

## Self-Check: PASSED

All 8 created files found. All 3 commits verified. All acceptance criteria met (14 tests total: 8 embed + 6 CLI).
