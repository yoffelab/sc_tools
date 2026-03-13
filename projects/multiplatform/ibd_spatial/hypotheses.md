---
type: hypotheses
project: ibd_spatial
platform: multiplatform
---

# ibd_spatial Hypotheses

Scientific questions for the IBD Spatial multiplatform project (51 samples, CosMx 1k/6k + Xenium MT/5K/colon, patient-matched).

## H1: Cross-platform integration quality

**Status:** #hypothesis/investigating
**Question:** Can patient-matched CosMx and Xenium data be integrated while preserving biological signal (CD vs UC, Inflamed vs Non-inflamed)?
**Approach:** Per-patient cross-platform mixing score (target >0.8 for Milestone 0). Gene intersection strategy with reduced latent dims (n_latent=10, n_hidden=64).

## H2: Platform batch effects vs biological signal

**Status:** #hypothesis/investigating
**Question:** Are platform-specific batch effects separable from disease biology in patient-matched samples?
**Approach:** Block ID (patient) is the biological anchor, NOT a batch to correct. Evaluate across 4 milestones (M0: technical baseline, M1: 16 matched patients, M2: high-plex, M3: full integration).

## H3: IBD biology recovery

**Status:** #hypothesis/blocked
**Question:** Can we recover known CD vs UC transcriptional differences from integrated cross-platform data?
**Blocker:** 6 of 7 panels have zero-filled RDS files (S3 sync failure). Only CosMx_6k usable. Contact Saha lab (jip2007) to re-sync from S3 or request AWS credentials.
