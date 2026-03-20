# Research Summary: sc_tools Agent-Native CLI

**Domain:** Agent-native CLI layer for computational biology pipelines
**Researched:** 2026-03-20
**Overall confidence:** HIGH

## Executive Summary

The agent-native CLI layer for sc_tools should use **Typer** (v0.24.1) as the CLI framework, **Pydantic** (v2.7+) for structured output schemas, and **custom JSON sidecar files** for provenance tracking. This stack adds only two new dependencies (typer, which brings rich/shellingham) while leveraging Pydantic already present in the ecosystem.

The key architectural insight is that "agent-native" means inverting the traditional CLI pattern: JSON is the default output (to stdout), and human-readable output is opt-in (via `--human`, rendered by Rich to stderr). Every command returns a Pydantic-validated `CLIResult` envelope so agents can reliably check status, extract data, and find output artifacts.

For provenance, the research strongly supports starting with lightweight JSON sidecars (`.provenance.json` alongside each output file) rather than formal standards like W3C PROV or RO-Crate. This aligns with the project's explicit decision to let the domain model emerge from usage before committing to a schema. The migration path is clean: custom sidecars now, W3C PROV-JSON export later, RO-Crate packaging for publication eventually.

Self-description (`sct schema`, `sct list-commands`, `sct describe <cmd>`) is a key differentiator. No widely-adopted standard exists for CLI self-description (OpenAPI is HTTP-only), but the pattern of exporting a JSON schema document -- combining Typer's command tree introspection with Pydantic's `model_json_schema()` -- gives agents a complete contract without parsing help text.

## Key Findings

**Stack:** Typer 0.24 + Pydantic 2.7 + Rich 13 + custom JSON provenance sidecars. Two new dependencies.
**Architecture:** JSON-to-stdout, Rich-to-stderr pattern. CLIResult Pydantic envelope for all commands. Provenance as file sidecars.
**Critical pitfall:** Mixing human and machine output on stdout will break agent parsing. Strict stdout=JSON, stderr=human separation is essential.

## Implications for Roadmap

Based on research, suggested phase structure:

1. **Foundation** - CLIResult envelope, output machinery, `sct` entry point
   - Addresses: Structured output format, dual-mode output (JSON/human)
   - Avoids: Building commands before the output contract is stable

2. **Core Commands** - Wrap existing sc_tools operations (`sct qc`, `sct preprocess`, `sct validate`)
   - Addresses: Table-stakes CLI commands, migration from ad-hoc scripts
   - Avoids: Starting with complex commands (ingest has multi-modality churn risk)

3. **Self-Discovery** - `sct schema`, `sct list-commands`, `sct describe`
   - Addresses: Agent self-discovery, machine-readable command catalog
   - Avoids: Building before commands exist (schema export needs commands to introspect)

4. **Provenance** - JSON sidecar writer, `sct provenance show/trace`
   - Addresses: File-based lineage tracking
   - Avoids: Premature DB migration

5. **Extended Commands** - `sct ingest`, `sct benchmark`, `sct mudata`
   - Addresses: Full pipeline coverage
   - Avoids: API churn (ingest modalities, benchmark formats still evolving)

**Phase ordering rationale:**
- Foundation must come first (all commands depend on CLIResult)
- Core commands validate the output contract before self-discovery exports it
- Provenance is additive (can be wired into existing commands incrementally)
- Extended commands depend on upstream sc_tools stabilization (benchmark pre-computed support)

**Research flags for phases:**
- Phase 5 (Extended Commands): Likely needs deeper research -- benchmark pre-computed format and MuData assembly patterns are still evolving
- Phase 1-3: Standard patterns, unlikely to need additional research

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack (CLI framework) | HIGH | Typer is dominant, well-maintained, version verified on PyPI |
| Stack (output format) | HIGH | JSON+Pydantic is industry standard (AWS CLI, gh, kubectl all use JSON) |
| Stack (provenance) | HIGH | Custom sidecars match project constraints; migration path to standards is clear |
| Stack (self-description) | MEDIUM | No standard exists; custom JSON schema export is pragmatic but unproven at scale |
| Features | HIGH | Derived directly from PROJECT.md requirements + industry CLI patterns |
| Architecture | HIGH | JSON-stdout/Rich-stderr pattern well-established in modern CLI tools |
| Pitfalls | HIGH | Grounded in codebase analysis (existing argparse scripts, MCP patterns) |

## Gaps to Address

- **Typer + Pydantic v2 integration patterns:** Typer's internal Pydantic usage may conflict with explicit Pydantic v2 models for output. Need to verify during implementation.
- **Shell completion for complex types:** Typer's completion may not handle AnnData-specific path patterns (`.h5ad` files). May need custom completers.
- **Provenance sidecar performance:** Writing checksums (sha256) for 25G h5ad files takes time. May need async or optional checksumming.
- **MCP-CLI symmetry:** Some operations exist in both MCP and CLI. Need a clear policy on which is source-of-truth (recommendation: CLI wraps sc_tools directly, MCP wraps CLI or sc_tools independently).

## Sources

- [Typer PyPI](https://pypi.org/project/typer/) -- v0.24.1, verified 2026-03-20
- [prov PyPI](https://pypi.org/project/prov/) -- v2.1.1, W3C PROV-DM
- [ro-crate-py GitHub](https://github.com/ResearchObject/ro-crate-py) -- RO-Crate 1.2
- [jc CLI tool](https://kellyjonbrazil.github.io/jc/) -- JSON conversion patterns
- [CLI framework comparison](https://dasroot.net/posts/2025/12/building-cli-tools-python-click-typer-argparse/)
- [RO-Crate provenance paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0309210)

---
*Research summary: 2026-03-20*
