# Plan: Visium HD QC Report Upgrade

**Scope:** Upgrade `pre_filter` (and `post_filter`) QC reports to support viable-region assessment —
the primary question being: "can we recover a usable tissue subregion from a failing sample?"

**Trigger:** Robin dataset reveals 3 samples (P3S3, P4S4, P6S6) with median counts 1-3 at 8µm
but 1.0–1.3M total counts and 13-14K genes detected. These are not failed samples — they are
samples where background bins outnumber tissue bins, collapsing the median. The current report
cannot distinguish this from true RNA degradation.

**Principle:** Do NOT save a filtered AnnData. All computations are report-only. The same
viable-region metrics must appear in the post-filter report so the user can see how they shift.

---

## Files to change

| File | Type of change |
|------|---------------|
| `sc_tools/qc/sample_qc.py` | Add `viable_threshold` param to `compute_sample_metrics` |
| `sc_tools/qc/plots.py` | Add `qc_viable_spots_spatial` and `qc_sample_viable_summary` |
| `sc_tools/qc/report.py` | Thread viable metrics into `generate_pre_filter_report` and `generate_post_filter_report` |
| `sc_tools/assets/pre_filter_qc_template.html` | Add viable-region table columns + section |
| `sc_tools/assets/post_filter_qc_template.html` | Same viable-region section (to show shift) |

---

## Step 1 — `compute_sample_metrics` additions

### Two-tier threshold design

The viable/tissue assessment uses two distinct thresholds with different purposes:

| Tier | Default | Name | Purpose |
|------|---------|------|---------|
| Tier 1 | counts ≥ 10, genes ≥ 5 | **tissue-present** | Confirms a bin is on tissue, not empty glass. Used for stroma detection (Step 1b). |
| Tier 2 | counts ≥ 100, genes ≥ 20 | **TPM-worthy** | Minimum for TPM normalization, signature scoring, cell typing at spot level. Used for recovery decisions. |

These are computed independently. A spot can be tissue-present (tier 1) but not TPM-worthy (tier 2).

**Rationale for 100/20 as TPM floor at 8µm:** A 55µm Visium spot (passing threshold ~500 counts) covers ~47x the area of an 8µm bin. Scaling: 500/47 ≈ 10 counts/bin at proportional density. 100 counts at 8µm implies a well-captured single-nucleus region. At 100 total counts / 20 genes, Poisson CV on the total is ~10% — acceptable for TPM. Individual gene estimates remain noisy but signature scoring aggregates over many spots.

Add parameters to `compute_sample_metrics`:
- `viable_threshold: dict | None = None` — Tier 1, default `{"min_counts": 10, "min_genes": 5}`
- `tpm_threshold: dict | None = None` — Tier 2, default `{"min_counts": 100, "min_genes": 20}`
- `soft_min_counts: int = 2` — lower bound for contextual (stromal) spots
- `neighborhood_k: int = 8` — number of spatial nearest neighbors for context computation
- `bin_area_um2: float | None = None` — area of one observation in µm². If None, inferred from modality (see table below). Pass 0 or False to disable area computation.

**Bin area by modality:**

| Modality | `bin_area_um2` | Notes |
|----------|---------------|-------|
| `visium_hd` | 64.0 | 8µm × 8µm bin |
| `visium_hd_cell` | None | Cell segmentation — variable area, skip |
| `visium` | 2376.0 | π × 27.5² µm (55µm diameter spot) |
| `xenium`, `cosmx` | None | Cell segmentation — variable area, skip |
| `imc` | None | Cell segmentation — variable area, skip |

### Step 1a — Tissue-present spots (Tier 1 viable, existing)

For each sample, apply an internal in-memory spot filter (do NOT modify `adata`):

```
viable_mask = (counts >= min_counts) & (n_genes >= min_genes)
```

New columns (rename from prior implementation to clarify tier):

| Column | Meaning |
|--------|---------|
| `n_viable_spots` | # spots passing the Tier 1 tissue-present filter |
| `pct_viable_spots` | `n_viable_spots / n_spots * 100` |
| `viable_total_counts_median` | Median counts among tissue-present spots |
| `viable_n_genes_median` | Median genes among tissue-present spots |
| `viable_pct_mt_median` | Median %MT among tissue-present spots (if MT computed) |
| `viable_pct_mt_gt20` | % tissue-present spots with >20% MT |

### Step 1a2 — TPM-worthy spots (Tier 2)

Apply the TPM threshold independently (not a subset of tissue-present — computed from all spots):

```
tpm_mask = (counts >= tpm_min_counts) & (n_genes >= tpm_min_genes)
```

New columns:

| Column | Meaning |
|--------|---------|
| `n_tpm_spots` | # spots passing the Tier 2 TPM-worthy filter |
| `pct_tpm_spots` | `n_tpm_spots / n_spots * 100` |
| `tpm_area_um2` | `n_tpm_spots * bin_area_um2` (None if bin_area_um2 is None) |
| `tpm_area_mm2` | `tpm_area_um2 / 1e6` (None if bin_area_um2 is None) |
| `tpm_n_genes_median` | Median genes among TPM-worthy spots |
| `tpm_pct_mt_median` | Median %MT among TPM-worthy spots (if MT computed) |
| `tpm_pct_mt_gt50` | % TPM-worthy spots with >50% MT (degradation in usable spots) |

The `tpm_area_mm2` column is the primary recovery metric: it directly answers "how much usable tissue area does this sample have?"

**Reference areas (Visium HD 8µm):**
- 500 TPM-worthy spots = 0.032 mm² (marginal)
- 1,000 spots = 0.064 mm² (minimum recommended for spatial analysis)
- 5,000 spots = 0.32 mm² (comfortable for most analyses)
- 10,000 spots = 0.64 mm² (good quality sample subregion)

### Step 1b — Contextual (stromal) spots (neighborhood-aware)

**Motivation:** Stroma, connective tissue, and ECM are biologically real but transcriptionally
quiet (~2-8 counts/bin at 8µm). A hard count threshold misclassifies them as background.
The discriminator is spatial context: a low-count bin surrounded by high-count tissue bins
is stroma. A low-count bin surrounded by other low-count bins in an empty region is background.

**Algorithm:** For each sample where `.obsm["spatial"]` is present:

1. Build k-NN graph from physical spatial coordinates (`obsm["spatial"]`), k=`neighborhood_k`
   — use `sklearn.neighbors.NearestNeighbors` with `algorithm="ball_tree"`, `metric="euclidean"`.
   Do NOT use squidpy (may not be available in all environments).
2. For each spot, compute `neighborhood_median_counts` = median of the `counts` values of its
   k nearest neighbors (excluding the spot itself).
3. Classify each spot:
   - **active**: counts ≥ `min_counts` (hard threshold, same as Step 1a)
   - **contextual**: counts ≥ `soft_min_counts` AND counts < `min_counts`
     AND `neighborhood_median_counts` ≥ `min_counts`
   - **background**: everything else

New columns:

| Column | Meaning |
|--------|---------|
| `n_contextual_spots` | # spots classified as contextual (stromal/low-signal tissue) |
| `pct_contextual_spots` | `n_contextual_spots / n_spots * 100` |
| `contextual_n_genes_median` | Median genes among contextual spots |
| `n_total_tissue_spots` | `n_viable_spots + n_contextual_spots` |
| `pct_total_tissue` | `n_total_tissue_spots / n_spots * 100` |

**Fallback:** If `.obsm["spatial"]` is absent (non-spatial modality like scRNA-seq), skip
Step 1b entirely — all contextual columns are `None`. The report gracefully omits the
contextual columns when they are all None.

**Performance note:** For large samples (100K+ spots), use a random subsample of 50K spots
to fit the k-NN index, then predict for all spots. This keeps runtime under 10 seconds.

Add `failure_mode` and `recovery_note` to `classify_samples`. Add parameters:
- `tpm_min_area_mm2: float = 0.05` — minimum TPM-worthy area to consider a sample recoverable
  (0.05 mm² ≈ 781 spots at 64 µm²/spot). Samples below this are classified as `insufficient_area`.
- `modality: str = "visium_hd"` — used to suppress MT-based rules for protein modalities (IMC).

**`failure_mode` rules (applied in order, first match wins):**

| Failure mode | Rule | Interpretation |
|-------------|------|----------------|
| `"complete_failure"` | `total_counts_sum < 0.1M` | Near-zero RNA recovery. No tissue. |
| `"rna_degradation"` | MT available AND `tpm_pct_mt_gt50 > 0.5` (>50% of TPM-worthy spots have >50% MT) | Even usable spots are predominantly mtRNA. Not recoverable. |
| `"insufficient_area"` | `tpm_area_mm2 < tpm_min_area_mm2` (and not above) | TPM-worthy spots exist but cover too little area for meaningful analysis. |
| `"sparse_tissue"` | `total_counts_mean_median_ratio > 5` AND `total_counts_sum >= 0.5M` AND `tpm_area_mm2 >= tpm_min_area_mm2` | Background dominates globally, but a recoverable tissue subregion exists. |
| `"low_transcription"` | fails QC threshold but none of above | Uniformly low signal. Manual review. |
| `""` | `qc_pass == True` | Passing sample. |

**Key design decisions:**
- MT threshold for `rna_degradation` is 50% (not 20%) in TPM-worthy spots. A spot with 100 counts and 40% MT still contributes 60 nuclear transcripts — analyzable. Only flag degradation when TPM-worthy spots are predominantly mitochondrial.
- `rna_degradation` only fires if MT data is present AND the TPM-worthy subpopulation itself is degraded. If MT not computed (e.g. no MT genes in panel), skip this rule entirely — never assign `rna_degradation` without evidence.
- `insufficient_area` is a soft failure — the sample has real RNA but not enough contiguous usable area. Default threshold 0.05 mm². This is configurable: a study needing large spatial context might raise it to 0.1 mm².
- For IMC (protein modality): skip `rna_degradation` rule entirely (no MT channel). `modality` param provides the explicit guard.

**`recovery_note` strings:**

| failure_mode | recovery_note |
|-------------|---------------|
| `complete_failure` | `"Near-zero transcript recovery. Not recoverable."` |
| `rna_degradation` | `f"RNA degradation in usable spots ({tpm_pct_mt_gt50*100:.0f}% of TPM-worthy spots >50% MT). Not recoverable by spatial filtering."` |
| `insufficient_area` | `f"TPM-worthy area too small: {tpm_area_mm2:.3f} mm² ({n_tpm_spots:,} spots). Below minimum {tpm_min_area_mm2} mm²."` |
| `sparse_tissue` | `f"Tissue subregion recoverable: {tpm_area_mm2:.3f} mm² ({n_tpm_spots:,} TPM-worthy spots). Consider spatial subsetting."` |
| `low_transcription` | `"Uniformly low transcription. Manual review recommended."` |
| `""` | `""` |

---

## Step 2 — New plots

### `qc_viable_spots_spatial(adata, metrics, sample_col, viable_threshold, soft_min_counts, neighborhood_k, ...)`

For each sample where spatial coordinates exist in `.obsm["spatial"]`:

**Three-class coloring (draw order: gray → amber → green, so tissue sits on top):**
- `#2ecc71` **green**: active tissue (counts ≥ min_counts)
- `#f39c12` **amber**: contextual/stromal (counts ≥ soft_min but < min_counts, neighborhood median ≥ min_counts)
- `#e0e0e0` **light gray**: background

Title: `f"{sample}\n{n_active:,} active + {n_contextual:,} stromal / {n_total:,} total"`
Subtitle: `f"MT% (active): {viable_pct_mt_median:.1f}%"` if MT available, else omit.
Add a small legend (3 colored patches: Active tissue / Stromal-contextual / Background).

Use `ax.scatter` with `s=max_point_size`, `linewidths=0`, `rasterized=True`.
No axis ticks or labels. Tight layout.
Return dict: `{sample_name: Figure}`.

Caller decides which samples to include. Include all by default.

### `qc_sample_viable_summary(metrics, classified, ...)`

Stacked bar chart per sample (log Y scale):
- Bottom stack: `n_viable_spots` (green `#2ecc71`)
- Middle stack: `n_contextual_spots` (amber `#f39c12`)
- Remaining: background (gray `#e0e0e0`)
- Total bar height = `n_spots`

Annotate `pct_total_tissue` (active + contextual) above each bar.
Sort samples by `pct_total_tissue` descending. Failed samples shown with hatched outline.
Legend: Active tissue / Stromal-contextual / Background.
X-tick labels: sample names, rotated 45°.
Title: `f"Tissue coverage per sample (active≥{min_counts} counts, stromal≥{soft_min} counts + neighbor context)"`

---

## Step 3 — Template changes

### Redundancies to remove (both pre and post templates)

- Remove "Failed Samples" bullet list below the table — duplicate of table `fail_reasons` column
- Remove `comparison_bar` (linear scale) — keep only `comparison_bar_log`
- Remove `violin_grouped` (linear scale) — keep only `violin_grouped_log`
- Remove `violin` (cohort-level pooled violin) — made redundant by `violin_grouped`
- In `qc_2x2`: drop linear `total_counts` panel; replace with a blank or `pct_viable_spots` bar

### Table additions (pre_filter_qc_template.html)

Add columns after existing columns, before "Fail reasons". Use `terms.observations_lower` for
observation labels (e.g., "cells" for xenium, "spots" for visium_hd):

```
| TPM-worthy {{ terms.observations_lower }} | % TPM | TPM area (mm²) | Tissue-present | % tissue | Stromal | Active %MT | %MT>50% TPM | Failure mode |
```

Column details:
- **TPM-worthy spots/cells**: `n_tpm_spots` — primary recovery metric
- **% TPM**: `pct_tpm_spots`
- **TPM area (mm²)**: `tpm_area_mm2` — shown only if `has_area` (spatial modalities with fixed bin size); format `"%.3f mm²"`
- **Tissue-present**: `n_viable_spots` (Tier 1)
- **% tissue**: `pct_viable_spots`
- **Stromal**: `n_contextual_spots` — shown only if `has_contextual`
- **Active %MT**: `tpm_pct_mt_median` — shown only if `has_mt`; format `"%.1f%%"`
- **%MT>50% TPM**: `tpm_pct_mt_gt50 * 100` — shown only if `has_mt`; format `"%.0f%%"`
- **Failure mode + recovery note**: combined column with `failure_mode — recovery_note`

`has_area` context variable: True if `tpm_area_mm2` is not None for any sample.

Show for all samples. Passing samples will show high TPM area, which provides positive reference.

### New section: "Viable Region Assessment"

Place this section between "Per-Sample Table" and "QC Distributions".

Contents:
1. `qc_sample_viable_summary` stacked bar chart (active + stromal + background per sample)
2. Per-sample spatial maps from `qc_viable_spots_spatial` in a 2-column grid, failed samples first.
   Caption: sample name, failure_mode, `n_active + n_stromal / n_total`, MT% if available.
3. Interpretive note:
   - "Green = active tissue (counts ≥ X, genes ≥ Y). Amber = stromal/contextual tissue
     (counts ≥ Z, neighborhood median ≥ X). Gray = background."
   - "High MT% in active spots indicates RNA degradation (not recoverable). Coherent
     green+amber region with low MT% indicates recoverable sample via spatial subsetting."

### Post-filter template

Add the same "Viable Region Assessment" section using the post-filter AnnData.
The metrics shift will reveal how effective the QC filter was:
- `pct_viable_spots` should increase post-filter (background removed)
- `viable_pct_mt_median` should be stable or decrease (background MT=0 removed)
- Samples that were borderline-fail may now show higher viable % → inform re-evaluation

---

## Step 4 — `generate_pre_filter_report` wiring

Add `viable_threshold` parameter (default `{"min_counts": 10, "min_genes": 5}`).

Steps:
1. Call `compute_sample_metrics` with `viable_threshold` to get enhanced metrics DataFrame
2. Call `classify_samples` as before (viable columns flow through automatically)
3. Generate `qc_viable_spots_spatial` → store as `plots["viable_spatial"]` dict
4. Generate `qc_sample_viable_summary` → store as `plots["viable_summary"]`
5. Pass `viable_threshold` into template context for display in interpretive note

Do the same for `generate_post_filter_report`.

**Do NOT write a filtered checkpoint. All filtered subsets exist only in memory during plot
generation.**

---

## Step 5 — MCP / `generate_qc_report` MCP tool

`mcp/tools/generate_qc_report` wraps `generate_pre_filter_report`. Add `viable_min_counts` and
`viable_min_genes` kwargs that flow through. Default values match Step 1 defaults.

---

## What NOT to change

- `apply_qc_filter` — do not touch, this is the actual filter that writes AnnData
- `filter_spots` — do not touch
- `classify_samples` pass/fail logic — do not change thresholds; add `failure_mode` as a
  supplemental column only, not a new pass/fail criterion
- The viable region is purely diagnostic/reporting — it informs human decisions, it does not
  automatically override QC pass/fail

---

## Implementation order

1. `compute_sample_metrics` — Step 1a (active viable) + Step 1b (contextual/stromal via k-NN)
2. `classify_samples` — add `failure_mode` and `recovery_note` columns; use `n_total_tissue_spots` in `sparse_tissue` rule
3. `report_utils.py` — `build_metrics_table_rows` passes all new columns through
4. `plots.py` — update `qc_viable_spots_spatial` (3-class colors) + update `qc_sample_viable_summary` (stacked bars)
5. `pre_filter_qc_template.html` — updated table columns + section + redundancy removal
6. `post_filter_qc_template.html` — same viable section + redundancy removal
7. `report.py` — add `soft_min_counts` and `neighborhood_k` params alongside existing viable params
8. MCP tool kwargs

---

## Success criteria

Success is defined as a **reviewer score ≥ 4.0 / 5** across all review dimensions below.
Two independent reviewers must both score ≥ 4. Anything below triggers a REVISE cycle.

---

### 1. Correct failure mode classification (must-pass, non-negotiable)

The report must distinguish the 5 Robin failed samples into exactly 2 categories:

**True failures — not recoverable (2 samples):**

| Sample | Expected `failure_mode` | Reason visible in report |
|--------|------------------------|--------------------------|
| `Pat5_Samp3_PRT_1_iE10_S6` | `complete_failure` | `total_counts_sum ≈ 0.02M`, `n_viable_spots ≈ 0`, spatial map shows no coherent region. Complete transcript recovery failure — specimen likely detached or handling failure. |
| `PT05-4_TUM` | `rna_degradation` | `viable_pct_mt_gt20` elevated (Plan.md documents 18.6% mean MT, up to 100%). Even tissue-covered spots show dying-cell signature. Not recoverable by spatial filtering. |

**Sparse tissue — potentially recoverable (3 samples):**

| Sample | Expected `failure_mode` | Key metrics visible in report |
|--------|------------------------|-------------------------------|
| `P3S3_1_D1_iG5_S13` | `sparse_tissue` OR `insufficient_area` | mean/median ≈ 6x, 1.24M total counts, 14K genes. `tpm_area_mm2` computed and shown. Spatial map shows coherent subregion. Recovery note includes area in mm². |
| `P4S4_1_A1_iB6_S16` | `sparse_tissue` OR `insufficient_area` | mean/median ≈ 10x, 1.02M total counts. Cell-seg version has 40K cells — spatial map must show clustered viable region. |
| `P6S6_2_A1_iB5_S8` | `sparse_tissue` OR `insufficient_area` | mean/median ≈ 7x, 1.26M total counts. Same pattern as P3S3. |

Note: which of `sparse_tissue` vs `insufficient_area` fires depends on actual `tpm_area_mm2`
computed from the data. Both are acceptable — the key requirement is that `tpm_area_mm2` is
computed, displayed, and the recovery note references the actual area value.

**`PT05-4_TUM` reclassification:** Now classified as `rna_degradation` only if
`tpm_pct_mt_gt50 > 0.5` in TPM-worthy spots. If that condition is not met (MT data unavailable
or most TPM-worthy spots are below 50% MT), it may be classified as `sparse_tissue` or
`insufficient_area` — which is acceptable, since it would then show a small but real
`tpm_area_mm2` that the user can evaluate. The threshold prevents false `rna_degradation`
labels when evidence is insufficient.

**The report must state WHY each sample is classified as it is**, including the `tpm_area_mm2`
value and `n_tpm_spots` count in the recovery note.

---

### 2. Spatial plot quality (reviewer-scored, weight: 30%)

- Spatial viable-spot maps for all 5 failed samples must be present and readable
- Maps use 3 colors: green (active), amber (stromal/contextual), gray (background)
- P3S3/P4S4/P6S6: coherent green+amber tissue subregion visible, spatially clustered
- Pat5: near-empty map (no green, no amber, only sparse gray)
- PT05-4_TUM: green spots present but report flags high MT% in active region → rna_degradation
- Stacked summary bar chart shows stromal (amber) contribution alongside active (green) per sample
- Caption includes: sample name, failure_mode, n_active + n_stromal / n_total, MT% if available

**Score 5:** All 5 maps clear, 3-class coloring correct, stacked summary bar present, captions complete
**Score 4:** All 5 maps present, minor caption or scale issues, stacked bar present
**Score 3:** Maps present but contextual/stromal class missing (only 2 colors) or stacked bar absent
**Score < 3:** Missing maps, wrong coloring, or active/stromal conflated → REVISE

---

### 3. Modality specificity (reviewer-scored, weight: 25%)

The report must adapt to modality. Reviewer verifies by running against a non-transcriptomics
modality (IMC or visium_hd_cell) or by code inspection of `get_modality_terms` usage:

- **Transcriptomics (visium_hd, visium_hd_cell, xenium, cosmx):**
  - MT% columns shown in viable region table
  - `viable_pct_mt_median` and `viable_pct_mt_gt20` used for `rna_degradation` classification
  - Labels read "genes", "spots/cells", "counts"

- **Protein/cytometry (IMC):**
  - MT% columns hidden (protein panels have no mitochondrial markers)
  - `failure_mode` classification uses only count/feature-based rules (no MT threshold)
  - Labels read "proteins", "cells", "intensity"
  - Viable region plot still shown (spatial intensity coverage still meaningful)

- **Cell segmentation (visium_hd_cell, xenium):**
  - Observation label is "cells" not "spots"
  - Viable threshold defaults should be lower (cells vs bins have different count scales)

**Score 5:** All three modality branches tested or code-verified correct
**Score 4:** Transcriptomics correct; IMC branch verified by code inspection
**Score 3:** Only transcriptomics tested, IMC untested but code path exists
**Score < 3:** MT% shown for IMC, or terminology is wrong for any tested modality → REVISE

---

### 4. Redundancy removal (reviewer-scored, weight: 20%)

Verify the following are absent from the rendered HTML:

- `[ ]` "Failed Samples" bullet list below the metrics table
- `[ ]` Linear-scale `comparison_bar` plot (only log10 version remains)
- `[ ]` Linear-scale `violin_grouped` plot (only log10 version remains)
- `[ ]` Pooled cohort-level `violin` plot (per-sample version makes it redundant)
- `[ ]` Linear `total_counts` panel in 2x2 grid (log1p version retained)

**Score 5:** All 5 removed, no regressions in remaining plots
**Score 4:** 4/5 removed
**Score 3:** 3/5 removed
**Score < 3:** Fewer than 3 removed, or a retained plot is broken → REVISE

---

### 5. Asset and rendering integrity (reviewer-scored, weight: 25%)

- All HTML template files in `sc_tools/assets/` render without Jinja2 errors
- `pre_filter_qc_template.html` and `post_filter_qc_template.html` both have the viable section
- `generate_pre_filter_report` and `generate_post_filter_report` both produce valid HTML
- No AnnData files written anywhere during report generation
- `compute_sample_metrics` returns viable columns when `viable_threshold` is set
- `classify_samples` returns `failure_mode` column
- Existing tests pass (`pytest sc_tools/tests/ -x -q`)
- New functions have docstrings and are exported in `__all__`

**Score 5:** All checks pass, new tests added for viable metric computation
**Score 4:** All checks pass, no new tests but existing pass
**Score 3:** HTML renders but one of pre/post is missing viable section
**Score < 3:** Jinja2 errors, AnnData written, or existing tests broken → REVISE

---

### Overall pass threshold

| Reviewer | Minimum score |
|----------|--------------|
| Reviewer 1 (code quality + correctness) | ≥ 4.0 |
| Reviewer 2 (report content + scientific value) | ≥ 4.0 |

If either reviewer scores < 4.0, or if Section 1 (failure mode classification) fails for ANY of
the 5 named samples, the implementation does not pass regardless of aggregate score.
