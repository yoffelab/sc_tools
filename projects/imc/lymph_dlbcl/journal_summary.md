# Journal Summary: lymph_dlbcl (Two-Panel IMC DLBCL)

Condensed summary. Full entries in `Journal.md`.

## Project scope

Two-panel Hyperion IMC on DLBCL (diffuse large B-cell lymphoma). Immune (T2) and stromal (S2) panels kept separate. Manuscript reproduction: 5 main + 8 supplementary = 13 figures total.

## Current status (2026-03-04)

- **Pipeline implemented:** 25 scripts + Snakefile + config.yaml for full manuscript reproduction
- **Approach:** Hybrid — reuse 47 Seurat-converted h5ad objects (on cayuga); consolidate into sc_tools checkpoints
- **Remote data:** `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2` (cayuga)
- **48 h5ad files** in `results/seurat_converted/` (immune T1/T2, stromal S1/S2, merged, spatial)
- **Inventory complete:** 1,668 remote files catalogued with annotated descriptors
- **Next:** Run download + phase0-1 on cayuga via `snakemake --cores 8 phase0 phase1`

## Key decisions

- **Panels:** T2 (immune, 1.63M cells) + S2 (stromal, 1.55M cells) as primary — most complete with DLC_code and labels
- **LME classes (5):** Cold (k-means 0,7,9), CD206 Enriched (1,8), Cytotoxic (2,10), Stroma (3,4), T cell Regulated (5,6) — from DLBCL_case_clustering.ipynb k=10 k-means
- **Checkpoint naming:** `adata.{immune|stromal}.{raw.p1|annotated.p2|celltyped.p4}.h5ad`
- **Spatial:** `adata.stromal.spatial.communities.h5ad` from SO_k30_community_cluster.h5ad
- **ML (Fig 5):** RandomForest + RandomizedSearchCV, seed=42, 5-fold CV
- **Stats:** BH FDR always, lifelines for survival, significance bars per skills.md

## Pipeline phases

| Phase | Scripts | Status |
|-------|---------|--------|
| 0 | download_clinical_metadata.sh, validate_h5ad_objects.py | Scripts ready |
| 1 | build_panel_adata.py, attach_clinical_metadata.py, build_spatial_adata.py | Scripts ready |
| 2-3 | validate_celltypes.py, build_lme_classes.py | Scripts ready |
| 4 | fig1-5 + supp_fig1-8 (13 scripts) | Scripts ready |
| 5 | Snakefile + config.yaml | Done |
| 6 | validate_figures.py | Script ready |

## Risks

1. S2 objects lack DLC_code — need cell ID CSVs for sample mapping
2. Spatial coords may differ between Seurat/AnnData
3. ML reproducibility depends on random seed
4. R-to-Python porting — validate numerically, not pixel-perfect
