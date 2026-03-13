# Clinical Data Schema

Registry table `subjects` stores cross-project de-identified patient/subject records. Each subject can link to multiple projects via `subject_project_links` and to multiple samples via the `samples` table.

---

## Subject Table Fields

| Column | Type | Required | Default | Description |
|--------|------|----------|---------|-------------|
| `id` | Integer | auto | auto-increment | Internal primary key |
| `subject_id` | String | yes | — | Unique de-identified identifier (cohort prefix + opaque ID, e.g. `DRAKE_001`) |
| `organism` | String | no | `human` | Species: `human`, `mouse`, `rat` |
| `sex` | String | no | — | `M`, `F`, or `unknown` |
| `age_at_collection` | Float | no | — | Age in years at time of tissue collection (integer or float, never exact DOB) |
| `diagnosis` | String | no | — | Free-text diagnosis (e.g. `prostate adenocarcinoma`, `DLBCL`) |
| `diagnosis_code` | String | no | — | Coded diagnosis: ICD-10 (e.g. `C61`) or OncoTree code |
| `disease_stage` | String | no | — | Clinical or pathological stage; staging system varies by disease (e.g. `pT3a`, `Stage IIIA`) |
| `treatment_status` | String | no | — | `treatment_naive`, `on_treatment`, or `post_treatment` |
| `tissue_of_origin` | String | no | — | Primary tissue site (e.g. `prostate`, `lung`, `pancreas`) |
| `vital_status` | String | no | — | `alive`, `deceased`, or `unknown` |
| `survival_days` | Float | no | — | Overall survival in days from diagnosis to last follow-up or death |
| `cause_of_death` | String | no | — | Free text if `vital_status = deceased` (e.g. `disease_progression`, `unrelated`) |
| `clinical_metadata_json` | Text | no | — | JSON overflow for cohort-specific fields (see below) |
| `created_at` | String | auto | UTC timestamp | Record creation time |

---

## clinical_metadata_json Overflow Pattern

Fields that do not fit the standard columns are stored as a JSON text blob in `clinical_metadata_json`. Parse with `json.loads()` in Python.

**When to use overflow vs. a standard column:**

- Standard columns cover fields shared across most cohorts (sex, age, diagnosis, staging, survival).
- Overflow handles cohort-specific fields that apply to one disease or study design.

**Example overflow content:**

```json
{
  "gleason_score": "4+3",
  "psa_at_diagnosis": 12.5,
  "smoking_status": "former",
  "pack_years": 30,
  "tumor_markers": {"CA19-9": 450, "CEA": 8.2},
  "recurrence_sites": ["liver", "peritoneum"],
  "prior_therapies": ["FOLFIRINOX", "gemcitabine"]
}
```

**Rules for overflow content:**
- Must be valid JSON (dict at top level)
- Must not contain PHI (see De-identification Rules below)
- Keys should use `snake_case`
- Numeric values stored as numbers, not strings

---

## Per-Cohort Column Mappings

| Cohort | Disease | Prefix | n | Sex | Key Overflow Fields | Notes |
|--------|---------|--------|---|-----|---------------------|-------|
| DRAKE | Prostate cancer | `DRAKE_` | 456 | All M | `gleason_score`, `psa_at_diagnosis`, `erg_status` | Radical prostatectomy cohort |
| ROS | Prostate cancer | `ROS_` | 135 | All M | `gleason_score`, `psa_at_diagnosis` | Radiation oncology series |
| GGO-IMC | Lung (GGO) | `GGO_IMC_` | 146 | Mixed | `smoking_status`, `pack_years`, `histology_subtype` | Ground-glass opacity cohort, IMC imaging |
| PLM | PDAC | `PLM_` | 38 | Mixed | `ca19_9`, `recurrence_sites`, `prior_therapies` | Survival recorded in months in source data; convert to days (multiply by 30.44) before storing in `survival_days` |

**Mapping source columns to schema:**

Each cohort has its own raw clinical spreadsheet with non-standard column names. When ingesting, map source columns to the standard `subjects` fields:

```python
# Example: DRAKE cohort mapping
column_map = {
    "Patient_ID": "subject_id",       # prefix with DRAKE_
    "Gender": "sex",                   # recode Male -> M
    "Age": "age_at_collection",
    "Path_Stage": "disease_stage",
    "Vital_Status": "vital_status",    # recode Living -> alive, Deceased -> deceased
    "OS_months": "survival_days",      # multiply by 30.44
    "Gleason": "gleason_score",        # -> clinical_metadata_json
}
```

---

## De-identification Rules

| Rule | Detail |
|------|--------|
| No PHI | No names, dates of birth, MRNs, addresses, phone numbers, or social security numbers |
| Age as number | Store age as integer or float in `age_at_collection`; never store exact date of birth |
| Opaque subject IDs | Subject IDs use cohort prefix + numeric or alphanumeric code (e.g. `DRAKE_042`); no initials or name fragments |
| Collection dates | Use relative dates (days from enrollment) or omit; never store exact calendar dates that could identify a subject |
| Overflow JSON | `clinical_metadata_json` must not contain any identifiable information; review before ingestion |
| Minimum cell size | If a subgroup has fewer than 5 subjects, consider suppressing or aggregating demographic breakdowns in public reports |

---

## Related Tables

- **`samples`** — physical specimens; links subject to project-scoped data (`tissue`, `fixation_method`, `batch`)
- **`subject_project_links`** — many-to-many join; `role` field: `enrolled` (default), `reference`, or `control`
- **`bio_data`** — typed biological data objects linked to samples and subjects

See [[registry]] for the full registry schema and BioData hierarchy.
