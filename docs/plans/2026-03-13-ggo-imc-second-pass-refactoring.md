---
title: ggo-imc Second-Pass Refactoring
project: ggo_imc
date: 2026-03-13
status: approved
---

# ggo-imc: Second-Pass Refactoring Plan (Final)

**Repo:** `projects/imc/ggo-imc/` | **Branch from:** `main`
**Context:** First-pass (`refactor/script-deduplication`) merged. This plan closes correctness, reproducibility, and maintainability gaps by aligning with scanpy/squidpy/sc_tools conventions.

---

## Hard Constraints (every phase, every agent)

1. **Figure outputs bit-identical** — no PDF filename, axis value, or color changes.
2. **Scientific logic untouched** — marker lists, statistical tests, score dicts preserved verbatim.
3. **Scripts remain runnable as `python scripts/foo.py`**
4. **Snakefile not touched**
5. **NEVER re-run clustering, embedding, or cell-typing steps.** Do not call or add:
   - `sc.tl.leiden`, `sc.tl.louvain`
   - `sc.tl.umap`, `sc.tl.tsne`
   - `sc.pp.neighbors`
   - Any clustering, neighbor graph, or embedding computation on data that has already been processed.

   These steps produce results that can change across runs (software versions, random seeds, hardware). Re-running them would invalidate the existing cell type annotations and require a full re-annotation cycle. All pipeline scripts load pre-computed `.h5ad` files with embeddings and annotations already stored — **these must not be overwritten or recomputed.** The `if os.path.exists(checkpoint):` guards in `myeloid_analysis.py` are an example of the correct pattern — preserve all such guards.

---

## Phases Summary

| # | Scope | Risk | Files |
|---|-------|------|-------|
| 1 | Remove dead subset code in `t_cell_analysis.py` | Low | 1 |
| 2 | Replace `print()` progress with `sc.logging.info` | Low | 3 |
| 3 | Consolidate duplicated `CYTOKINE` constant into `utils.py` | Low | 3 |
| 4 | Decision: `patient_group.py` lines 121-122 left untouched | None | 0 |
| 5 | Smoke-test harness (`tests/test_imports.py`) | Low | 1 new |
| 6 | Pytest unit tests for `cond_prob` + `load_config` | Medium | 2 new |
| 7 | New plot utility abstractions | **Deferred** | — |

Phases 1-4 unconditional. Phases 5-7 require team sign-off.

---

## Cross-Cutting Verification

- **Phases 1-5:** Import smoke-test only — `python -c 'from scripts.utils import load_config, load_panels, save_figure, ensure_dir'`. Do not run pytest (directory does not exist until Phase 6).
- **Phase 6+:** `pytest projects/imc/ggo-imc/tests/ -v`
- **Never** run a full pipeline script as a verification step — doing so risks writing new `.h5ad` files over pre-computed checkpoints.

---

## Phase 1 — Dead Code Removal (`t_cell_analysis.py`)

Lines 21-26 build a PANEL_H lymphocyte subset immediately overwritten by `sc.read()` at lines 28-31.
Lines 52-58 repeat the same pattern for PANEL_G, overwritten at lines 60-63.

**Lines to delete:**
- Lines 21-26: `h_lymphocyte_markers` definition + dead `h_lymphocyte_idx` / `h_lymphocytes` assignments
- Lines 52-58: `g_lymphocyte_markers` definition + dead `g_lymphocyte_idx` / `g_lymphocytes` assignments

**Lines to keep (LIVE):** 28-31 and 60-63 (the `sc.read()` calls that load pre-computed lymphocyte AnnData).

**`myeloid_analysis.py` excluded from Phase 1.** Its subset construction at lines 30-33 and 67-70 feeds an `if os.path.exists(...)` branch that runs UMAP/neighbor computation only when the cached `.h5ad` is absent. This is live code protecting against recomputation — do not touch it.

**Commit:** `fix(t_cell_analysis): remove dead AnnData subset code overwritten by sc.read`

---

## Phase 2 — Progress Logging (`sc.logging.info`)

Replace informational progress `print()` with `sc.logging.info`. Scientific output prints are not touched (see Phase 4).

**Targets:**
- `epithelial_characterization.py` lines 23, 68: reading progress → `sc.logging.info`
- `utils.py` line 93 inside `load_panels`: `print(f"Reading {path}...")` → `sc.logging.info`
- `ue_analysis.py`: no standalone progress prints (tqdm handles iteration). No changes.

Note: `sc.logging.info` writes to stderr. No Snakemake rule in this project uses a `log:` directive, so stdout/stderr distinction is inconsequential.

**Commit:** `refactor(logging): replace print() with sc.logging.info in progress messages`

---

## Phase 3 — Constants Consolidation (`utils.py`)

Only constants with 2+ identical occurrences move to `utils.py`.

| Constant | Source | Lines | Why |
|----------|--------|-------|-----|
| `CYTOKINE` | `t_cell_analysis.py` | 33 | Identical list in `myeloid_analysis.py:55`. Two occurrences. |

`FUNCTIONAL_MARKERS` (`ue_analysis.py`, single occurrence) and `PROCESSES` (`patient_group.py`, single occurrence) stay inline per the 2+ rule.

**Files touched:** `utils.py`, `t_cell_analysis.py`, `myeloid_analysis.py`

**Commit:** `refactor(constants): consolidate repeated biology constants into utils.py`

---

## Phase 4 — Decision: Leave `patient_group.py` Lines 121-122 Unchanged

`patient_group.py` lines 121-122 print computed DataFrames representing `P(histology | radiology)` — scientific results, not progress messages. Source lines are left entirely unchanged. This decision is recorded here only; no edit to the source file is required.

---

## Phase 5 — Smoke-Test Harness

Create `tests/test_imports.py` with one bare `import` test per pipeline script. No data access required — `if __name__ == '__main__'` guards all data loading.

**Commit:** `test: add import smoke-test harness for pipeline scripts`

---

## Phase 6 — Pytest Unit Tests (Team Gate)

Add unit tests for pure functions only:
- `patient_group.cond_prob` — small synthetic DataFrame; verify pivot shape and row sums to 100
- `utils.load_config` — fixture YAML; verify dict keys

Do not test any function that loads real AnnData or calls `sc.tl.*` / `sc.pp.*`.

**Gate:** Team sign-off required before starting.

**Commit:** `test: add unit tests for cond_prob and utils functions`

---

## Phase 7 — Deferred Indefinitely

- `annotated_regplot`: single occurrence in `t_cell_analysis.py`. Fails 2+ rule.
- `label_stacked_bar`: appears in `epithelial_characterization.py` and `patient_group.py` but patterns are NOT identical (label-centering math differs). A shared abstraction would require conditional logic obscuring intent.

No justified work. Phase 7 deferred until a genuine 2+ identical pattern is identified.

---

## Sequencing

```
Phase 1 → Phase 2 → Phase 3 → Phase 5 → Phase 6 (team gate)
Phase 4: decision only, no source change, any time.
Phase 7: deferred.
```

---

## Critical Files

| File | Phases | Notes |
|------|--------|-------|
| `scripts/t_cell_analysis.py` | 1 | Delete lines 21-26 and 52-58; keep 28-31 and 60-63 |
| `scripts/utils.py` | 2, 3 | Add `CYTOKINE` constant; change logging in `load_panels` |
| `scripts/myeloid_analysis.py` | 3 | Import `CYTOKINE` from utils; EXCLUDED from Phase 1 |
| `scripts/patient_group.py` | 4 | Lines 121-122 untouched |
| `scripts/epithelial_characterization.py` | 2 | Logging at lines 23, 68 |
