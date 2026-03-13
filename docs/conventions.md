---
type: conventions
tier: tool-design
tags:
  - convention
---

# Conventions

Coding, statistical, and naming standards for all sc_tools projects. Canonical source: `skills.md` at repo root.

Full standards (statistics, figures, lint): [[skills]]

## Naming (wiki-specific)

- **Signature scores:** `obsm['signature_score']` / `obsm['signature_score_z']`; not in `obs` by default
- **Colors:** `adata.uns[f'{name}_colors']`; never overwrite if length matches
- **Checkpoint names:** See [[checkpoints.gen]] for required metadata per phase
- **Paths in scripts:** use `$PROJECT_DIR` or `$(PROJECT)` — never hardcode repo-root paths
