---
name: k-dense-comm
tier: 3
description: "[Tier 3 -- Specialized] Scientific communication reference: IMRaD drafting, journal figure specs, grant structures, peer review responses, citation formatting, and reporting guidelines. Use when the science-writer agent needs concrete templates, checklists, or formatting rules for publication or grant work."
allowed-tools: Read Write Edit Bash Glob Grep WebSearch WebFetch
---

# Scientific Communication

## Core Principle

Write in flowing prose, not bullets. Outline first (research-lookup), then convert to paragraphs.

## IMRaD Section Templates

### Introduction
1. Establish importance (broad to specific, 2-3 paragraphs)
2. Review literature and identify gap
3. State hypothesis or research question
4. Describe approach and novelty (1 paragraph)

### Methods -- Reproducibility Checklist
- [ ] Software with version numbers (e.g., scanpy 1.9.6, scvi-tools 1.0.4)
- [ ] All parameter values for non-default settings
- [ ] Statistical test name, assumptions checked, correction method
- [ ] Sample sizes per group and how determined
- [ ] Ethical approvals and consent statements
- [ ] Data/code availability with accession numbers or repository URLs

### Results
- Report findings objectively; no interpretation
- Every claim backed by: test statistic, effect size, adjusted p-value, N
- Reference figures/tables by number; integrate, do not duplicate

### Discussion
- Paragraph 1: principal findings restated
- Compare with literature (agreements and conflicts)
- Mechanistic speculation clearly labeled as such
- Limitations paragraph (do not bury)
- Future directions (1 paragraph max)

### Abstract
- Write LAST, as one flowing paragraph (unless journal requires structured)
- 150-250 words: background (1-2 sentences), methods, results, conclusion

## Journal Figure Requirements

| Spec | Nature | Cell | Science |
|------|--------|------|---------|
| Single col width | 89 mm | 85 mm | 55 mm |
| Double col width | 183 mm | 178 mm | 175 mm |
| Min font size | 5 pt | 6 pt | 6 pt |
| Font family | Helvetica/Arial | Arial | Helvetica/Arial |
| Panel labels | **a**, **b** lowercase bold | **A**, **B** uppercase bold | **A**, **B** uppercase bold |
| Min DPI (raster) | 300 | 300 | 300 |
| Vector formats | PDF, EPS | PDF, EPS | PDF, EPS |
| Color model | RGB | RGB | RGB |

### Figure Legend Structure
1. Title: one-sentence summary of finding
2. Panel descriptions: (**a**) what is shown, N, test, P-value
3. Define error bars (SD, SEM, 95% CI), scale bars, color palette
4. Statistical details: test name, correction method, significance thresholds

### Colorblind-Safe Defaults
- Categorical: Okabe-Ito `['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']`
- Sequential: viridis, plasma, cividis (never jet/rainbow)
- Diverging: RdBu_r, PuOr (never red-green)

## Citation Styles

| Style | In-text | Key format rule |
|-------|---------|----------------|
| Nature | Superscript^1^ | Abbreviated journal, bold volume, year in parens |
| APA 7 | (Smith et al., 2024) | Author-date, ampersand before last author |
| Vancouver | [1] | Numbered brackets, no periods after initials |
| Cell | (Smith et al., 2024) | "and" before last author, no issue number |

BibTeX keys: `FirstAuthor2024keyword`. Protect caps with `{}`. Use `--` for page ranges.

## Grant Writing Structures

### NIH Specific Aims Page (1 page)
```
Paragraph 1: Hook + gap + long-term goal + central hypothesis
Aim 1: [Action verb] ... Working hypothesis: ...
Aim 2: [Action verb] ... Working hypothesis: ...
Aim 3: [Action verb] ... Working hypothesis: ...
Payoff paragraph: Expected outcomes + impact on field
```
- Aims must be independent (failure of one does not sink the project)
- Each aim: rationale, hypothesis, brief approach, expected outcome

### NIH Research Strategy (12 pages R01)
- **Significance**: important problem, rigor of prior research, how knowledge will advance
- **Innovation**: novel concepts, approaches, methodologies
- **Approach**: preliminary data, research design, methods, alternatives, timeline

### NSF Project Description (15 pages)
- Broader Impacts integrated throughout (not an appendix)
- Prior NSF support section with results and publications
- Intellectual Merit: potential to advance knowledge

### NIH Review Scoring (1-9; 1 = exceptional)
Significance | Investigator(s) | Innovation | Approach | Environment

### NSF Review Criteria (equally weighted)
Intellectual Merit | Broader Impacts

## Peer Review Response Template

```
## Reviewer #N

> [Quoted reviewer comment]

**Response:** We thank the reviewer for this observation. [Address point with evidence.]
We have revised [section/page/line]. The updated text now reads: "[new text]".
[If new analysis:] As shown in new Supplementary Figure X, ...
```
- Address every point; never dismiss
- For disagreements: provide evidence, remain respectful
- Track changes summary at end: list all modified sections

## Reporting Guidelines (Key Items)

**STROBE** (observational): study design in title, eligibility criteria, confounders, flow diagram.
**PRISMA** (systematic reviews): registered protocol, search strategy per DB, flow diagram, risk-of-bias.
**Computational biology**: code repo + version tag, container spec, raw data accession, random seeds.

## Hypothesis Formulation

Observation -> Mechanism -> Prediction -> Falsification -> Experimental design.
Generate 3-5 competing hypotheses; rank by parsimony and testability.

## Spatial Transcriptomics Methods Boilerplate

**Visium HD**: "Tissue sections were placed on Visium HD slides (10x Genomics) and processed per manufacturer protocol (CG000450). Libraries were sequenced on Illumina NovaSeq 6000. Reads were aligned to [GRCh38/GRCm39] using Space Ranger v[X.Y.Z]. QC and analysis used scanpy v[X.Y.Z] and squidpy v[X.Y.Z]."

**Deconvolution**: "Cell type deconvolution was performed using [Cell2location/Tangram/RCTD] with a reference scRNA-seq dataset from [source, accession], filtered to [N] cell types (min [N] cells per type). Results stored in `obsm['cell_type_proportions']`."

## Common Writing Errors

| Error | Fix |
|-------|-----|
| "significant" without p-value | Always report: P = 0.003, FDR-corrected |
| Causal language from correlation | "associated with" not "caused by" |
| Gene symbol not italicized | *YTHDF2* (gene), YTHDF2 (protein) |
| Bar chart for distributions | Use violin + strip plot |
| "data not shown" | Show it or remove the claim |
| Passive voice throughout | Mix active/passive; prefer active for clarity |
