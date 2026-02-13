# Journal Summary: lymph_dlbcl (Two-Panel IMC DLBCL)

Condensed summary of `Journal.md` for this project. Full entries are in `Journal.md`.

## Project scope

- Two-panel Hyperion IMC (DLBCL): **immune** and **stromal** panels, kept separate (e.g. separate h5ad per panel). Data is on remote path; priority is to **reuse existing processed data** from the previous analyst; process raw only when needed.

## Recent phase (2025-02-12)

- **Inventory:** Remote file listing (CSV, h5ad, rds, h5); script `build_remote_inventory.py` produces `metadata/remote_file_inventory.csv`; notebook scan adds data_type, analysis_stage, needs_download. **Annotated** version with natural-language descriptors: `scripts/annotate_inventory_descriptors.py` writes `metadata/remote_file_inventory_annotated.csv` (detailed_descriptor filled). Original script output kept unchanged for traceability.
- **Required files:** Normalized expression (cells as observations) per panel is primary; raw expression per panel when available. Identify in inventory → `metadata/files_to_download.csv` → download to `data/downloaded/` → `metadata/download_manifest.csv`. Build AnnData with normalized (and raw in layers when present).
- **RDS:** `analyze_rds_notebook_usage.py` produced `metadata/rds_notebook_usage.csv` and `metadata/notebook_rds_summary.csv` for traceability of 47 RDS files to notebooks.

## Key conventions

- Panels: **immune**, **stromal** (separate checkpoints).
- Remote (read-only): `/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2`.
- Use annotated inventory for download decisions; keep script output and annotated file separate.
