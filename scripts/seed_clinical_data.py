"""Seed the sc_tools registry with real clinical data from XLSX files.

Parses clinical spreadsheets for DRAKE, ROS, GGO-IMC, and PLM cohorts,
then populates Subject and Sample tables with standardized fields.

Usage:
    python scripts/seed_clinical_data.py [--dry-run] [--replace-synthetic]
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent

# ── XLSX source definitions ──────────────────────────────────────────────────

DRAKE_ROS_XLSX = REPO_ROOT / "projects" / "imc" / "aapc" / "metadata" / "Hyperion ROS-DRAKE.xlsx"
GGO_IMC_XLSX = (
    REPO_ROOT
    / "projects"
    / "imc"
    / "ggo-imc"
    / "metadata"
    / "De-identified J&J Clinical Annotations_UpdatedAAH031523 (1).xlsx"
)
PLM_XLSX = (
    REPO_ROOT
    / "projects"
    / "imc"
    / "pancreas_liver_meta"
    / "metadata"
    / "15-015_TMAset_clinicopath_data.xlsx"
)

# ── Project name mappings (for linking subjects to registry projects) ────────

COHORT_PROJECT_MAP: dict[str, str] = {
    "DRAKE": "aapc",
    "ROS": "aapc",
    "GGO_IMC": "ggo_human",
    "PLM": "pancreas_liver_meta",
}

# Synthetic subject ID prefixes seeded by seed_biodata.py
SYNTHETIC_PREFIXES = [
    "DLBCL_",
    "GGO_V_",
    "GGO_IMC_",
    "ROBIN_",
    "IBD_",
    "AAPC_",
    "CLONEVO_",
    "PLM_",
    "BCG_",
    "HL_",
]


# ── Helpers ──────────────────────────────────────────────────────────────────


def _nan_to_none(val: object) -> object:
    """Convert pandas NaN/NaT to None for clean JSON and DB storage."""
    if pd.isna(val):
        return None
    return val


def standardize_vital_status(val: object) -> str:
    """Normalize vital status to alive / deceased / unknown."""
    if val is None or pd.isna(val):
        return "unknown"
    s = str(val).strip().lower()
    if s in ("alive", "ned", "ned/alive", "living"):
        return "alive"
    if s in ("dead", "dod", "deceased", "duk", "died", "died of disease"):
        return "deceased"
    return "unknown"


def standardize_sex(val: object) -> str:
    """Normalize sex to M / F / unknown."""
    if val is None or pd.isna(val):
        return "unknown"
    s = str(val).strip().lower()
    if s in ("m", "male", "1"):
        return "M"
    if s in ("f", "female", "2"):
        return "F"
    return "unknown"


def _build_json_overflow(row: pd.Series, keys: list[str]) -> str | None:
    """Build a JSON string from selected row keys, skipping NaN values."""
    overflow: dict[str, object] = {}
    for k in keys:
        if k in row.index:
            v = _nan_to_none(row[k])
            if v is not None:
                overflow[k] = v
    return json.dumps(overflow) if overflow else None


# ── Cohort parsers ───────────────────────────────────────────────────────────


def parse_drake(xlsx_path: Path) -> tuple[list[dict], list[dict]]:
    """Parse the DRAKE sheet. Returns (subjects, samples)."""
    df = pd.read_excel(xlsx_path, sheet_name="DRAKE")
    subjects: list[dict] = []
    samples: list[dict] = []

    overflow_keys = [
        "race",
        "psa_at_dx",
        "gleason_score",
        "grade_group",
        "surgical_t_stage",
        "n_stage",
        "recurrence",
        "time_to_recurrence",
        "treatment_after_recurrence",
        "adjuvant_treatment",
    ]

    # Build a mapping from actual columns to overflow keys (case-insensitive)
    col_map: dict[str, str] = {}
    for col in df.columns:
        col_lower = col.strip().lower().replace(" ", "_")
        for ok in overflow_keys:
            if ok == col_lower or ok in col_lower:
                col_map[col] = ok

    for _, row in df.iterrows():
        raw_id = _nan_to_none(row.get("id"))
        if raw_id is None:
            continue

        subject_id = f"DRAKE_{raw_id}"

        # Collect overflow from original column names mapped to clean keys
        overflow: dict[str, object] = {}
        for orig_col, clean_key in col_map.items():
            v = _nan_to_none(row.get(orig_col))
            if v is not None:
                overflow[clean_key] = v

        vital = standardize_vital_status(row.get("current status"))
        cod = _nan_to_none(row.get("death cause"))

        subjects.append(
            {
                "subject_id": subject_id,
                "organism": "human",
                "sex": "M",
                "diagnosis": "prostate_cancer",
                "tissue_of_origin": "prostate",
                "age_at_collection": _nan_to_none(row.get("Age at Surgery")),
                "disease_stage": _nan_to_none(row.get("Clinical T Stage")),
                "vital_status": vital,
                "cause_of_death": cod,
                "clinical_metadata_json": json.dumps(overflow) if overflow else None,
            }
        )

        # Samples from TMA grid
        tma = _nan_to_none(row.get("TMA#"))
        gridnum = _nan_to_none(row.get("tma_gridnum"))
        if tma is not None:
            sample_id = f"DRAKE_{raw_id}_T{tma}"
            if gridnum is not None:
                sample_id = f"{sample_id}_{gridnum}"
            samples.append(
                {
                    "sample_id": sample_id,
                    "subject_id": subject_id,
                    "tissue": "prostate",
                    "fixation_method": "FFPE",
                    "sample_type": "TMA",
                }
            )

    logger.info("DRAKE: parsed %d subjects, %d samples", len(subjects), len(samples))
    return subjects, samples


def parse_ros(xlsx_path: Path) -> tuple[list[dict], list[dict]]:
    """Parse the ROS sheet. Returns (subjects, samples)."""
    df = pd.read_excel(xlsx_path, sheet_name="ROS")
    subjects: list[dict] = []
    samples: list[dict] = []

    overflow_keys = [
        "race",
        "psa_at_dx",
        "gleason_score",
        "grade_group",
        "surgical_t_stage",
        "n_stage",
        "recurrence",
        "time_to_recurrence",
        "treatment_after_recurrence",
        "adjuvant_treatment",
    ]

    col_map: dict[str, str] = {}
    for col in df.columns:
        col_lower = col.strip().lower().replace(" ", "_")
        for ok in overflow_keys:
            if ok == col_lower or ok in col_lower:
                col_map[col] = ok

    for _, row in df.iterrows():
        raw_id = _nan_to_none(row.get("SUBJECT_ID"))
        if raw_id is None:
            continue

        subject_id = f"ROS_{raw_id}"

        overflow: dict[str, object] = {}
        for orig_col, clean_key in col_map.items():
            v = _nan_to_none(row.get(orig_col))
            if v is not None:
                overflow[clean_key] = v

        vital = standardize_vital_status(row.get("Vital Status"))
        cod = _nan_to_none(row.get("death cause"))

        subjects.append(
            {
                "subject_id": subject_id,
                "organism": "human",
                "sex": "M",
                "diagnosis": "prostate_cancer",
                "tissue_of_origin": "prostate",
                "age_at_collection": _nan_to_none(row.get("Age at Surgery")),
                "disease_stage": _nan_to_none(row.get("Clinical T Stage")),
                "vital_status": vital,
                "cause_of_death": cod,
                "clinical_metadata_json": json.dumps(overflow) if overflow else None,
            }
        )

        # TMA samples: tumor and benign cores
        tumor_tma = _nan_to_none(row.get("Tumor_TMA#"))
        benign_tma = _nan_to_none(row.get("Benign_TMA#"))
        if tumor_tma is not None:
            samples.append(
                {
                    "sample_id": f"ROS_{raw_id}_tumor",
                    "subject_id": subject_id,
                    "tissue": "prostate",
                    "tissue_region": "tumor",
                    "fixation_method": "FFPE",
                    "sample_type": "TMA",
                }
            )
        if benign_tma is not None:
            samples.append(
                {
                    "sample_id": f"ROS_{raw_id}_benign",
                    "subject_id": subject_id,
                    "tissue": "prostate",
                    "tissue_region": "benign",
                    "fixation_method": "FFPE",
                    "sample_type": "TMA",
                }
            )

    logger.info("ROS: parsed %d subjects, %d samples", len(subjects), len(samples))
    return subjects, samples


def parse_ggo_imc(xlsx_path: Path) -> tuple[list[dict], list[dict]]:
    """Parse the GGO-IMC clinical sheet. Returns (subjects, samples)."""
    df = pd.read_excel(xlsx_path, sheet_name="Sheet1")
    subjects: list[dict] = []
    samples: list[dict] = []

    overflow_keys = [
        "race",
        "smoking_status",
        "packyrs",
        "radiology",
        "solid_comp",
        "tumor_size",
        "molecular",
        "histology",
        "invasion",
        "margins",
        "lymph_nodes",
    ]

    col_map: dict[str, str] = {}
    for col in df.columns:
        col_lower = col.strip().lower().replace(" ", "_")
        for ok in overflow_keys:
            if ok == col_lower or ok in col_lower:
                col_map[col] = ok

    for row_idx, (_, row) in enumerate(df.iterrows()):
        subject_id = f"GGO_IMC_{row_idx}"

        overflow: dict[str, object] = {}
        for orig_col, clean_key in col_map.items():
            v = _nan_to_none(row.get(orig_col))
            if v is not None:
                overflow[clean_key] = v

        sex = standardize_sex(row.get("Gender"))

        subjects.append(
            {
                "subject_id": subject_id,
                "organism": "human",
                "sex": sex,
                "diagnosis": "GGO",
                "tissue_of_origin": "lung",
                "age_at_collection": _nan_to_none(row.get("Age")),
                "disease_stage": _nan_to_none(row.get("PATH STAGE")),
                "vital_status": "unknown",
                "clinical_metadata_json": json.dumps(overflow) if overflow else None,
            }
        )

        # GGO-IMC does not have explicit TMA/sample columns;
        # one subject maps to one sample (the ROI).
        samples.append(
            {
                "sample_id": f"GGO_IMC_{row_idx}_ROI",
                "subject_id": subject_id,
                "tissue": "lung",
                "fixation_method": "FFPE",
                "sample_type": "resection",
            }
        )

    logger.info("GGO-IMC: parsed %d subjects, %d samples", len(subjects), len(samples))
    return subjects, samples


def parse_plm(xlsx_path: Path) -> tuple[list[dict], list[dict]]:
    """Parse the PLM (pancreas-liver metastasis) clinical sheet. Returns (subjects, samples)."""
    df = pd.read_excel(xlsx_path, sheet_name="Human-CZ_coded2")
    subjects: list[dict] = []
    samples: list[dict] = []

    overflow_keys = [
        "differentiation",
        "t_stage",
        "n_stage",
        "lvi",
        "pni",
        "surgical_margin",
        "bmi",
        "diabetes",
        "neoadjuvant",
        "adjuvant",
        "tumor_size",
        "tumor_location",
        "ca19_9",
    ]

    col_map: dict[str, str] = {}
    for col in df.columns:
        col_lower = col.strip().lower().replace(" ", "_")
        for ok in overflow_keys:
            if ok == col_lower or ok in col_lower:
                col_map[col] = ok

    for _, row in df.iterrows():
        raw_id = _nan_to_none(row.get("ID"))
        if raw_id is None:
            continue

        subject_id = f"PLM_{raw_id}"

        overflow: dict[str, object] = {}
        for orig_col, clean_key in col_map.items():
            v = _nan_to_none(row.get(orig_col))
            if v is not None:
                overflow[clean_key] = v

        sex = standardize_sex(row.get("Sex"))

        # Convert OS from months to days
        os_months = _nan_to_none(row.get("OS"))
        survival_days = None
        if os_months is not None:
            try:
                survival_days = float(os_months) * 30.44
            except (TypeError, ValueError):
                survival_days = None

        vital = standardize_vital_status(row.get("Death_Status"))

        subjects.append(
            {
                "subject_id": subject_id,
                "organism": "human",
                "sex": sex,
                "diagnosis": "PDAC",
                "tissue_of_origin": "pancreas",
                "age_at_collection": _nan_to_none(row.get("Age")),
                "disease_stage": _nan_to_none(row.get("AJCC_8th_stage")),
                "survival_days": survival_days,
                "vital_status": vital,
                "clinical_metadata_json": json.dumps(overflow) if overflow else None,
            }
        )

        # TMA samples
        tma_val = _nan_to_none(row.get("TMA_1_or_2"))
        if tma_val is not None:
            sample_id = f"PLM_{raw_id}_TMA{tma_val}"
            samples.append(
                {
                    "sample_id": sample_id,
                    "subject_id": subject_id,
                    "tissue": "pancreas",
                    "fixation_method": "FFPE",
                    "sample_type": "TMA",
                }
            )

    logger.info("PLM: parsed %d subjects, %d samples", len(subjects), len(samples))
    return subjects, samples


# ── Registry operations ──────────────────────────────────────────────────────


def remove_synthetic_subjects(reg: object) -> int:
    """Delete subjects whose IDs match synthetic prefixes from seed_biodata.py.

    Returns the number of deleted subjects.
    """
    from sqlalchemy import or_

    deleted = 0
    with reg._session() as sess:
        filters = [reg._Subject.subject_id.like(f"{prefix}%") for prefix in SYNTHETIC_PREFIXES]
        synthetic = sess.query(reg._Subject).filter(or_(*filters)).all()
        for subj in synthetic:
            sess.delete(subj)
            deleted += 1
        sess.commit()
    logger.info("Deleted %d synthetic subjects", deleted)
    return deleted


def insert_subjects(
    reg: object,
    subjects: list[dict],
    cohort_name: str,
    *,
    dry_run: bool = False,
) -> int:
    """Insert subjects into the registry. Returns count of successfully added subjects."""
    count = 0
    project_name = COHORT_PROJECT_MAP.get(cohort_name)

    for rec in subjects:
        sid = rec["subject_id"]
        if dry_run:
            logger.info("[DRY RUN] Would add subject: %s", sid)
            count += 1
            continue
        try:
            # Separate standard kwargs from subject_id and organism
            kwargs = {k: v for k, v in rec.items() if k not in ("subject_id", "organism")}
            reg.add_subject(rec["subject_id"], organism=rec.get("organism", "human"), **kwargs)
            count += 1
        except Exception:
            logger.exception("Failed to add subject %s", sid)
            continue

        # Link to project if it exists in the registry
        if project_name:
            try:
                reg.link_subject_to_project(sid, project_name)
            except Exception:
                logger.debug("Could not link %s to project %s", sid, project_name)

    return count


def insert_samples(
    reg: object,
    samples: list[dict],
    cohort_name: str,
    *,
    dry_run: bool = False,
) -> int:
    """Insert samples into the registry. Returns count of successfully added samples."""
    count = 0
    project_name = COHORT_PROJECT_MAP.get(cohort_name)

    for rec in samples:
        sample_id = rec["sample_id"]
        if dry_run:
            logger.info("[DRY RUN] Would add sample: %s", sample_id)
            count += 1
            continue
        try:
            kwargs = {k: v for k, v in rec.items() if k not in ("sample_id", "subject_id")}
            reg.add_sample(
                sample_id,
                subject_id=rec.get("subject_id"),
                project_name=project_name,
                **kwargs,
            )
            count += 1
        except Exception:
            logger.exception("Failed to add sample %s", sample_id)

    return count


# ── Main ─────────────────────────────────────────────────────────────────────


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Seed registry with real clinical data from XLSX files."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be inserted without writing to the database.",
    )
    parser.add_argument(
        "--replace-synthetic",
        action="store_true",
        help="Delete synthetic subjects (from seed_biodata.py) before inserting real data.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    from sc_tools.registry import Registry

    reg = Registry()

    # Optionally remove synthetic placeholder subjects
    if args.replace_synthetic and not args.dry_run:
        print("\n=== Removing synthetic subjects ===")
        n_removed = remove_synthetic_subjects(reg)
        print(f"  Removed {n_removed} synthetic subjects (and cascaded samples/links)")

    # Define cohort parsers: (cohort_name, xlsx_path, parser_function)
    cohort_defs: list[tuple[str, Path, object]] = [
        ("DRAKE", DRAKE_ROS_XLSX, parse_drake),
        ("ROS", DRAKE_ROS_XLSX, parse_ros),
        ("GGO_IMC", GGO_IMC_XLSX, parse_ggo_imc),
        ("PLM", PLM_XLSX, parse_plm),
    ]

    summary: dict[str, dict[str, int]] = {}

    for cohort_name, xlsx_path, parser_fn in cohort_defs:
        print(f"\n=== {cohort_name} ===")

        if not xlsx_path.exists():
            print(f"  SKIP: XLSX not found at {xlsx_path}")
            logger.warning("XLSX not found for %s: %s", cohort_name, xlsx_path)
            summary[cohort_name] = {"subjects": 0, "samples": 0}
            continue

        try:
            subjects, samples = parser_fn(xlsx_path)
        except Exception:
            logger.exception("Failed to parse %s from %s", cohort_name, xlsx_path)
            summary[cohort_name] = {"subjects": 0, "samples": 0}
            continue

        n_subj = insert_subjects(reg, subjects, cohort_name, dry_run=args.dry_run)
        n_samp = insert_samples(reg, samples, cohort_name, dry_run=args.dry_run)
        summary[cohort_name] = {"subjects": n_subj, "samples": n_samp}

        print(f"  Subjects: {n_subj}, Samples: {n_samp}")

    # Print final summary
    print("\n=== Summary ===")
    total_subj = 0
    total_samp = 0
    for cohort_name, counts in summary.items():
        prefix = "[DRY RUN] " if args.dry_run else ""
        print(
            f"  {prefix}{cohort_name}: {counts['subjects']} subjects, {counts['samples']} samples"
        )
        total_subj += counts["subjects"]
        total_samp += counts["samples"]

    print(f"\n  Total: {total_subj} subjects, {total_samp} samples")

    if not args.dry_run:
        status = reg.status()
        print("\n=== Registry Totals ===")
        print(f"  Projects:  {status['n_projects']}")
        print(f"  Subjects:  {status['n_subjects']}")
        print(f"  Samples:   {status['n_samples']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
