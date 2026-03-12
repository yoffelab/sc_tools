# Project Architecture: Two-Panel IMC DLBCL (lymph_dlbcl)

This document describes the **project-specific** directory layout, data sources, and data flow for the DLBCL two-panel IMC project. Repository-wide conventions are in the `docs/Architecture.md`.

---

## 1. Directory Layout

```text
projects/imc/lymph_dlbcl/
├── data/              # Local copy of raw or staged data (when used)
│   └── downloaded/   # Required processed/normalized/filtered files pulled from remote (Step 2)
├── figures/           # QC, manuscript, and exploratory figures
├── metadata/          # Project metadata and remote inventory
│   ├── remote_file_inventory.csv           # Script output (Step 2a); detailed_descriptor empty
│   ├── remote_file_inventory_annotated.csv # Annotated inventory with detailed_descriptor (Step 2b); use for download decisions
│   ├── remote_csv_listing.txt     # Raw find output from remote (reproducible input)
│   ├── files_to_download_rds.csv  # List of RDS paths to download (47 files)
│   ├── rds_notebook_usage.csv     # Per-RDS: which notebooks read/write it, inferred content (analyze_rds_notebook_usage.py)
│   ├── notebook_rds_summary.csv  # Per-notebook: which RDS it reads/writes (analyze_rds_notebook_usage.py)
│   ├── files_to_download.csv     # Step 2: list of remote paths to download (required for AnnData)
│   ├── download_manifest.csv     # After download: remote_path, local_path, date (traceability)
│   ├── sample_metadata.csv        # metadata_attach: sample → clinical (when ready)
│   └── celltype_map.json         # celltype_manual: cluster_id → celltype (when ready)
├── notebooks/         # Existing DLBCL notebooks (reference + extraction)
│   └── dlbcl_notebooks/
│       └── DLBCLv2/   # Two-panel immune and stromal preprocessing and analyses
├── results/           # AnnData checkpoints (standard names per root Architecture.md)
├── outputs/           # Intermediate outputs (logs, temp)
├── scripts/           # Project-specific ingestion and pipeline scripts
├── Mission.md
├── Architecture.md   # This file
└── Journal.md
```

---

## 2. Data Sources and Paths

### Remote (reference; no bulk download required for planning)

| Purpose | Path |
|--------|------|
| **Remote project root** | `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2` |
| **Access** | SSH (e.g. `ssh cayuga`) then run `find` or list commands. |

All CSV and directory-structure discovery is done against this path. Notebooks from the prior analyst may reference other paths (e.g. `/athena/elementolab/scratch/...`, `/home/dym2001/...`); the **canonical** source for this project is the backup path above.

### Local

| Purpose | Path |
|--------|------|
| **Project root** | `projects/imc/lymph_dlbcl/` |
| **Notebooks (reference)** | `projects/imc/lymph_dlbcl/notebooks/dlbcl_notebooks/DLBCLv2/` |
| **Inventory and metadata** | `projects/imc/lymph_dlbcl/metadata/` |
| **Checkpoint AnnData** | `projects/imc/lymph_dlbcl/results/` |

---

## 3. Checkpoint Nomenclature (alignment with root Architecture.md)

Once the pipeline phases are implemented, use the standard checkpoint names:

| Phase | Standard path (under `results/`) | Description |
|-------|----------------------------------|-------------|
| `qc_filter` | `adata.filtered.h5ad` | Raw AnnData after ingestion and QC. |
| `metadata_attach` | `adata.annotated.h5ad` | With sample/clinical metadata attached. |
| `preprocess` | `adata.normalized.h5ad` | Filtered, normalized, clustered. |
| `scoring` | `adata.scored.h5ad` | Normalized + signature scores (and optional automated cell typing). |
| `celltype_manual` | `adata.celltyped.h5ad` | Manual cell type labels applied. |

For **two-panel** structure, keep **immune** and **stromal** panels **separate**. Use panel-prefixed checkpoints (e.g. `results/adata.immune.raw.p1.h5ad`, `results/adata.stromal.raw.p1.h5ad`, `results/adata.immune.normalized.p3.h5ad`, `results/adata.stromal.normalized.p3.h5ad`). Document any additional naming in this file and in Mission.md.

---

## 4. Data Flow (Planning and Execution)

**Priority:** Reuse existing processed data from the previous analyst when the inventory shows it is suitable. Process raw only when needed.

1. **Step 1 — Inventory:** Remote `find` for `*.csv` (and optionally `*.h5ad`, `*.rds`, `*.h5`); notebook scan for path references; build `metadata/remote_csv_inventory.csv` with descriptor, inferred_role, and processed/raw tag.
2. **Step 2 (priority) — Identify, download, reuse:** From inventory, identify **required** processed/normalized/filtered files; record in `metadata/files_to_download.csv`; **download** them to `data/downloaded/`; record in `metadata/download_manifest.csv`; build AnnData from downloaded files into `results/adata.immune.*.h5ad` and `results/adata.stromal.*.h5ad`.
3. **Step 3 (fallback) — Raw:** For panel(s) without suitable processed data, ingest raw and run `qc_filter`; write panel-specific `results/adata.{immune|stromal}.raw.p1.h5ad` and downstream checkpoints. Keep immune and stromal pipelines separate.

---

## 5. Snakemake Pipeline

The full pipeline is defined in `Snakefile` + `config.yaml`. Run with:

```bash
snakemake --cores 8 all          # Full pipeline
snakemake --cores 8 phase0       # Download + validate
snakemake --cores 8 phase1       # Build AnnData checkpoints
snakemake --cores 8 figures      # All 13 figures
bash scripts/run_on_cayuga.sh all  # Sync + run on cayuga
```

### DAG

```
download_clinical -> validate_h5ad
                         |
                   build_panel_adata (immune + stromal)
                         |
                   attach_clinical_metadata
                         |
              +----------+----------+
              |                     |
        validate_celltypes    build_lme_classes
              |                     |
        build_spatial         (P4 checkpoints)
              |                     |
              +----------+----------+
                         |
           fig1, fig2, fig3, fig4, fig5
           supp_fig1 .. supp_fig8
                         |
                    validate_figures
```

### Scripts inventory (25 total)

| Script | Phase | Purpose |
|--------|-------|---------|
| download_clinical_metadata.sh | 0 | Download clinical/metadata CSVs from cayuga backup |
| validate_h5ad_objects.py | 0 | Validate shape, columns, labels of key h5ad objects |
| build_panel_adata.py | 1 | Build immune/stromal P1 checkpoints from T2/S2 objects |
| attach_clinical_metadata.py | 1 | Join clinical data -> P2 checkpoints |
| build_spatial_adata.py | 1 | Build spatial community object |
| validate_celltypes.py | 2-3 | Marker expression heatmap per cell type |
| build_lme_classes.py | 3 | Assign 5 LME classes -> P4 checkpoints |
| fig1_single_cell_atlas.py | 4 | UMAP, heatmap, dotplot, proportions |
| fig2_lme_classes.py | 4 | LME heatmap, abundance, composition |
| fig3_clinical.py | 4 | KM curves, COO, mutations, Cox |
| fig4_spatial.py | 4 | Spatial communities, neighborhoods |
| fig5_ml_framework.py | 4 | RF classification, ROC, IHC validation |
| supp_fig1_qc_panels.py | 4 | Panel QC metrics |
| supp_fig2_bcell.py | 4 | B cell subclusters |
| supp_fig3_tcell_myeloid.py | 4 | T cell / myeloid |
| supp_fig4_vessel.py | 4 | Vessel analysis |
| supp_fig5_tme_sensitivity.py | 4 | TME clustering sensitivity |
| supp_fig6_mutations.py | 4 | Mutation landscape |
| supp_fig7_rna_protein.py | 4 | RNA-protein comparison |
| supp_fig8_extended_survival.py | 4 | Extended survival |
| validate_figures.py | 6 | Verify figures + numerical validation |
| run_on_cayuga.sh | - | Sync + run pipeline on cayuga |

---

## 6. References

- Root `Architecture.md`: checkpoint names, required metadata, phase details.
- Root `Mission.md`: toolkit phasing and entry points.
- This project `Mission.md`: roadmap and active tasks.
- `config.yaml`: Project parameters, input paths, LME definitions, ML config.
- `Snakefile`: Full pipeline DAG.
