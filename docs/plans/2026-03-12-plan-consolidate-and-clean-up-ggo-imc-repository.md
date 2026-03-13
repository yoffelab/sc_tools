---
status: active
created: 2026-03-12
project: ggo_imc
summary: "Clean up ggo-imc repo: archive exploratory scripts, update Snakefile, rewrite README as publication companion"
---
# Plan: Consolidate and Clean Up ggo-imc Repository

## Context

The ggo-imc repo (https://github.com/ElementoLab/ggo-imc) is the publication companion for:

> **"Spatial profiling of early-stage lung adenocarcinoma reveals patterns of immunomodulation and epithelial plasticity"**
> Kim, Ravichandran, Yoffe et al. — Cancer Cell submission (Dec 2025)

The study uses IMC to profile 2.24M cells across 122 early-stage LUAD specimens, identifying two progression trajectories (immune-inflamed vs fibrotic) and a diagnostic blind spot where 20.5% of fibrotic tumors are radiologically misclassified. The manuscript has 5 main figures + extended figures.

The repo has 41 scripts in `scripts/`, but only 10 Python + 1 R (`asd.R`, not yet in repo) are in the active Snakefile pipeline. The rest are exploratory, legacy, or duplicated. Goal: clean up file organization, remove dead code, and update documentation to clearly map code to manuscript figures. No script refactoring.

### Manuscript figure-to-script mapping

| Figure | Manuscript content | Script(s) |
|--------|-------------------|-----------|
| **Figure 1** | Cohort overview, cell type heatmap (1a-d), ROI PCA archetypes (1e) | `celltype_heatmap_info.py`, `roi_pca_plot.py` |
| **Figure 2** | Immune dynamics: lymphocyte abundance across stages (2a-h), myeloid/macrophage polarization (2i-m) | `celltype_differential_abundance.py` (lymphoid myeloid), `t_cell_analysis.py`, `myeloid_analysis.py` |
| **Figure 3** | Stromal expansion, epithelial remodeling, EMP (3a-g) | `celltype_differential_abundance.py` (stromal epithelial), `epithelial_characterization.py` |
| **Figure 4** | UTAG microenvironments, TLS, tumor-stroma interface, cell-cell interactions (4a-f) | `ue_analysis.py` |
| **Figure 5** | Patient risk groups via hierarchical clustering, fibrotic trajectory, diagnostic gap (5a-e) | `roi_pca_plot_group.py`, `patient_group.py`, `asd.R` (R, not yet in repo) |

---

## Phase 1: Archive dead code and backups

**Files to modify:** `scripts/` directory

**Delete backup scripts** (all exist in git history):
- `myeloid_analysis_backup.py`
- `roi_pca_plot_backup.py`
- `ue_analysis_backup.py`
- `t_cell_analysis_backp.py`

**Move 24 exploratory scripts to `scripts/exploratory/`:**
- `rapids_scanpy_funcs.py`, `cell_proportion_across_condition.py`, `patient_feature_selection.py`, `cell_phenotype_distance.py`, `yaml_to_csv.py`, `cell_phenotyping.py`, `cell_phenotyping.ipynb`, `celltype_interaction.py`, `celltype_comparison_functional.py`, `clinical_features.py`, `sample_similarity.py`, `myeloid_reanalysis.py`, `cell_interaction.py`, `spatial_plot.py`, `utag_ue.py`, `patient_cell_prop.py`, `wes_mutation_analysis.py`, `wgs_mutation_analysis.py`, `mutation.py`, `emt.py`, `emt_morphology.py`, `create_patient_density_matrix.py`, `patient_density_heatmap_annotation.Rmd`, `celltype_differential_smokers.py`

**Delete deprecated makefile** (Snakefile is the active pipeline).

**Remove junk:**
- `scripts/Missing Data Panels` (text file, not a directory)
- `scripts/__pycache__/`
- `scripts/.ipynb_checkpoints/`

**Note:** `scripts_backup/`, `previous_code/`, `figures_backup/`, `documentation/` are already gitignored and untracked — no git action needed. They stay as local-only directories.

**Resulting `scripts/` layout:**
```
scripts/
  # Pipeline (Figures 1-5)
  celltype_heatmap_info.py           # Figure 1: cell type heatmap
  roi_pca_plot.py                    # Figure 1: ROI PCA archetypes
  celltype_differential_abundance.py # Figures 2, 3: immune/stromal density
  t_cell_analysis.py                 # Figure 2: lymphocyte functional states
  myeloid_analysis.py                # Figure 2: myeloid/macrophage polarization
  epithelial_characterization.py     # Figure 3: epithelial phenotypes, EMP
  ue_analysis.py                     # Figure 4: UTAG microenvironments, TLS
  roi_pca_plot_group.py              # Figure 5: patient group PCA overlay
  patient_group.py                   # Figure 5: patient risk stratification
  # Utilities
  download_yaml.py                   # Download data from Zenodo/Box
  concat_anndata.py                  # Concatenate per-sample AnnData
  label_metadata_anndata.py          # Attach clinical metadata
  generate_manifests.py              # Generate batch manifests
  # Benchmarks
  benchmark_integration.py           # Integration method comparison
  benchmark_segmentation.py          # Segmentation benchmark
  # Archived
  exploratory/                       # 24 exploratory scripts (not in pipeline)
```

---

## Phase 2: Update Snakefile

**File:** `Snakefile`

- Add `clean` rule to remove `*.done` sentinel files
- Add comment on `patient` rule noting `asd.R` is not yet in the repo (the rule references it but the file is missing)

---

## Phase 3: Update README.md

**File:** `README.md`

Rewrite to serve as the publication companion landing page:
- Add manuscript title, abstract summary, and citation placeholder
- Keep Zenodo DOI badges
- Replace `make` instructions with Snakemake commands
- Update Python version (3.9 → 3.10+)
- Add repository structure section showing the clean layout
- Add figure-to-script mapping table referencing manuscript figure content
- Remove the HTML comment block at the bottom

---

## Phase 4: Cleanup odds and ends

**Files:** `.gitignore`, `CLAUDE.md`, `pyproject.toml`

- `.gitignore`: add `*.done`, `__pycache__/`, `.DS_Store`, `.ipynb_checkpoints/`
- `CLAUDE.md`: add `scripts/exploratory/` to key files table, note `asd.R` is not yet in repo, add manuscript reference
- `pyproject.toml`: already has `>=3.10` — no change needed

---

## Verification

After Phase 1, confirm all Snakefile-referenced scripts still exist:
```bash
# Check each script referenced in Snakefile exists
grep -oP 'scripts/\S+\.py' Snakefile | sort -u | while read f; do test -f "$f" && echo "OK: $f" || echo "MISSING: $f"; done
```

---

## Commit Strategy

One commit per phase:
1. `chore: archive exploratory scripts, remove backups and deprecated makefile`
2. `chore: add clean rule and asd.R note to Snakefile`
3. `docs: rewrite README as publication companion with figure-script mapping`
4. `chore: update gitignore and CLAUDE.md for new repo structure`

Then push to origin/main.
