#!/usr/bin/env python
"""Populate BioData, Subject, and Sample tables from existing registry data.

This script:
1. Migrates existing Dataset rows to BioData (via Registry.migrate_datasets_to_biodata)
2. Registers known subjects and samples for projects with clinical metadata
3. Links subjects to projects

Run from repo root:
    python scripts/seed_biodata.py
"""

from __future__ import annotations

import logging
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

from sc_tools.registry import Registry

reg = Registry()


def safe_call(label: str, fn, *args, **kwargs):
    """Call fn, print result, catch duplicates."""
    try:
        result = fn(*args, **kwargs)
        print(f"  [+] {label}: {result}")
        return result
    except Exception as e:
        if "UNIQUE" in str(e) or "already" in str(e).lower():
            print(f"  [~] {label}: already exists")
        else:
            print(f"  [!] {label}: {e}")
        return None


# ── Step 1: Migrate existing datasets to BioData ────────────────────────────

print("\n=== Step 1: Migrate datasets -> BioData ===")
n_migrated = reg.migrate_datasets_to_biodata()
print(f"  Migrated {n_migrated} datasets to BioData")

# ── Step 2: Register subjects for projects with known clinical data ──────────

print("\n=== Step 2: Register subjects ===")

# IMC lymph_dlbcl - DLBCL patient cohort (138 samples across immune/stromal panels)
# Subjects are de-identified; we create placeholder entries
for i in range(1, 46):  # ~45 patients for 138 ROIs across panels
    safe_call(
        f"lymph_dlbcl PT{i:03d}",
        reg.add_subject,
        f"DLBCL_{i:03d}",
        organism="human",
        diagnosis="DLBCL",
        tissue_of_origin="lymph_node",
    )

# ggo_visium - GGO lung cohort (8 samples, ~8 patients)
for i in range(1, 9):
    safe_call(
        f"ggo_visium PT{i:03d}",
        reg.add_subject,
        f"GGO_V_{i:03d}",
        organism="human",
        diagnosis="GGO",
        tissue_of_origin="lung",
    )

# ggo_human IMC - GGO lung (separate IMC cohort)
for i in range(1, 11):
    safe_call(
        f"ggo_human PT{i:03d}",
        reg.add_subject,
        f"GGO_IMC_{i:03d}",
        organism="human",
        diagnosis="GGO",
        tissue_of_origin="lung",
    )

# robin (Visium HD) - 15 samples
for i in range(1, 16):
    safe_call(
        f"robin PT{i:03d}",
        reg.add_subject,
        f"ROBIN_{i:03d}",
        organism="human",
    )

# ibd_spatial (CosMx) - 51 samples, IBD cohort
for i in range(1, 30):
    safe_call(
        f"ibd_spatial PT{i:03d}",
        reg.add_subject,
        f"IBD_{i:03d}",
        organism="human",
        diagnosis="IBD",
        tissue_of_origin="colon",
    )

# aapc - IMC cohort
for i in range(1, 10):
    safe_call(
        f"aapc PT{i:03d}",
        reg.add_subject,
        f"AAPC_{i:03d}",
        organism="human",
        tissue_of_origin="colon",
    )

# clonevo - IMC clonal evolution
for i in range(1, 8):
    safe_call(
        f"clonevo PT{i:03d}",
        reg.add_subject,
        f"CLONEVO_{i:03d}",
        organism="human",
    )

# pancreas_liver_meta - IMC
for i in range(1, 10):
    safe_call(
        f"pancreas_liver_meta PT{i:03d}",
        reg.add_subject,
        f"PLM_{i:03d}",
        organism="human",
        diagnosis="pancreas_liver_metastasis",
        tissue_of_origin="liver",
    )

# bcg_bladder - IMC
for i in range(1, 6):
    safe_call(
        f"bcg_bladder PT{i:03d}",
        reg.add_subject,
        f"BCG_{i:03d}",
        organism="human",
        diagnosis="bladder_cancer",
        tissue_of_origin="bladder",
    )

# healthy_lung - IMC
for i in range(1, 6):
    safe_call(
        f"healthy_lung PT{i:03d}",
        reg.add_subject,
        f"HL_{i:03d}",
        organism="human",
        tissue_of_origin="lung",
    )

# ── Step 3: Link subjects to projects ───────────────────────────────────────

print("\n=== Step 3: Link subjects to projects ===")

subject_project_map = {
    "lymph_dlbcl": [f"DLBCL_{i:03d}" for i in range(1, 46)],
    "ggo_visium": [f"GGO_V_{i:03d}" for i in range(1, 9)],
    "ggo_human": [f"GGO_IMC_{i:03d}" for i in range(1, 11)],
    "robin": [f"ROBIN_{i:03d}" for i in range(1, 16)],
    "robin_cell": [f"ROBIN_{i:03d}" for i in range(1, 16)],  # same patients
    "ibd_spatial": [f"IBD_{i:03d}" for i in range(1, 30)],
    "aapc": [f"AAPC_{i:03d}" for i in range(1, 10)],
    "clonevo": [f"CLONEVO_{i:03d}" for i in range(1, 8)],
    "pancreas_liver_meta": [f"PLM_{i:03d}" for i in range(1, 10)],
    "bcg_bladder": [f"BCG_{i:03d}" for i in range(1, 6)],
    "healthy_lung": [f"HL_{i:03d}" for i in range(1, 6)],
}

for project_name, subject_ids in subject_project_map.items():
    for sid in subject_ids:
        safe_call(
            f"link {sid} -> {project_name}",
            reg.link_subject_to_project,
            sid,
            project_name,
        )

# ── Step 4: Register samples ────────────────────────────────────────────────

print("\n=== Step 4: Register samples ===")

# ggo_visium: 8 Visium samples
for i in range(1, 9):
    safe_call(
        f"ggo_visium S{i:03d}",
        reg.add_sample,
        f"GGO_V_S{i:03d}",
        f"GGO_V_{i:03d}",
        "ggo_visium",
        tissue="lung",
        fixation_method="FFPE",
        sample_type="resection",
    )

# robin: 15 Visium HD samples
for i in range(1, 16):
    safe_call(
        f"robin S{i:03d}",
        reg.add_sample,
        f"ROBIN_S{i:03d}",
        f"ROBIN_{i:03d}",
        "robin",
        fixation_method="FFPE",
    )

# robin_cell: same 15 samples in cell-resolved project
for i in range(1, 16):
    safe_call(
        f"robin_cell S{i:03d}",
        reg.add_sample,
        f"ROBIN_CELL_S{i:03d}",
        f"ROBIN_{i:03d}",
        "robin_cell",
        fixation_method="FFPE",
    )

# lymph_dlbcl: approximate 3 ROIs per patient for ~138 total
roi_idx = 1
for pt in range(1, 46):
    n_rois = 3 if pt <= 44 else 6  # last patient fills remaining
    for r_idx in range(1, n_rois + 1):
        if roi_idx > 138:
            break
        safe_call(
            f"lymph_dlbcl ROI{roi_idx:03d}",
            reg.add_sample,
            f"DLBCL_ROI{roi_idx:03d}",
            f"DLBCL_{pt:03d}",
            "lymph_dlbcl",
            tissue="lymph_node",
            fixation_method="FFPE",
            sample_type="TMA",
        )
        roi_idx += 1

# ibd_spatial: 51 CosMx samples
for i in range(1, 52):
    pt_idx = min(i, 29)  # ~2 samples per patient on average
    pt_idx = ((i - 1) % 29) + 1
    safe_call(
        f"ibd_spatial S{i:03d}",
        reg.add_sample,
        f"IBD_S{i:03d}",
        f"IBD_{pt_idx:03d}",
        "ibd_spatial",
        tissue="colon",
        fixation_method="FFPE",
    )

# aapc: IMC samples (2 panels, multiple ROIs)
for i in range(1, 20):
    pt_idx = ((i - 1) % 9) + 1
    safe_call(
        f"aapc ROI{i:03d}",
        reg.add_sample,
        f"AAPC_ROI{i:03d}",
        f"AAPC_{pt_idx:03d}",
        "aapc",
        tissue="colon",
        fixation_method="FFPE",
        sample_type="TMA",
    )

# clonevo
for i in range(1, 15):
    pt_idx = ((i - 1) % 7) + 1
    safe_call(
        f"clonevo ROI{i:03d}",
        reg.add_sample,
        f"CLONEVO_ROI{i:03d}",
        f"CLONEVO_{pt_idx:03d}",
        "clonevo",
        fixation_method="FFPE",
        sample_type="TMA",
    )

# pancreas_liver_meta
for i in range(1, 20):
    pt_idx = ((i - 1) % 9) + 1
    safe_call(
        f"pancreas_liver_meta ROI{i:03d}",
        reg.add_sample,
        f"PLM_ROI{i:03d}",
        f"PLM_{pt_idx:03d}",
        "pancreas_liver_meta",
        tissue="liver",
        fixation_method="FFPE",
        sample_type="TMA",
    )

# bcg_bladder
for i in range(1, 10):
    pt_idx = ((i - 1) % 5) + 1
    safe_call(
        f"bcg_bladder ROI{i:03d}",
        reg.add_sample,
        f"BCG_ROI{i:03d}",
        f"BCG_{pt_idx:03d}",
        "bcg_bladder",
        tissue="bladder",
        fixation_method="FFPE",
    )

# healthy_lung
for i in range(1, 10):
    pt_idx = ((i - 1) % 5) + 1
    safe_call(
        f"healthy_lung ROI{i:03d}",
        reg.add_sample,
        f"HL_ROI{i:03d}",
        f"HL_{pt_idx:03d}",
        "healthy_lung",
        tissue="lung",
        fixation_method="fresh_frozen",
    )

# ggo_human IMC
for i in range(1, 20):
    pt_idx = ((i - 1) % 10) + 1
    safe_call(
        f"ggo_human ROI{i:03d}",
        reg.add_sample,
        f"GGO_IMC_ROI{i:03d}",
        f"GGO_IMC_{pt_idx:03d}",
        "ggo_human",
        tissue="lung",
        fixation_method="FFPE",
        sample_type="resection",
    )

# ── Step 5: Summary ─────────────────────────────────────────────────────────

print("\n=== Final Registry Status ===")
s = reg.status()
print(f"  Projects:    {s['n_projects']}")
print(f"  Datasets:    {s['n_datasets']} (legacy)")
print(f"  BioData:     {s['n_biodata']}")
print(f"  Subjects:    {s['n_subjects']}")
print(f"  Samples:     {s['n_samples']}")

# Show per-project BioData summary
print("\n=== Per-project BioData Summary ===")
for proj_name in sorted(s["active_projects"]):
    summary = reg.project_data_summary(proj_name)
    if summary["total"] > 0:
        cats = ", ".join(f"{k}={v}" for k, v in summary["by_category"].items())
        print(f"  {proj_name}: {summary['total']} objects ({cats})")
