---
status: in_progress
created: 2026-03-12
summary: Finish IMC segmentation benchmark on brb — retry 21 missing tasks, aggregate results, review per-dataset F1 stats
---

# IMC Segmentation Benchmark — Finish & Aggregate

## Context

Parallel SLURM array job 13974653 ran 112 ROIs x 3 Cellpose strategies = 336 runs on brb.
91/112 tasks produced CSVs. 21 failed (pip race condition + possible timeouts on large GGO ROIs).
Retry job 13992020 submitted for the 21 missing tasks.

## Datasets

| Dataset | ROIs | Path |
|---------|------|------|
| lung_covid | 58 | `/athena/elementolab/scratch/thb4002/imc/tiffs/lung-true-512/processed` |
| ggo_imc | 54 | `/athena/elementolab/scratch/liy2010/IMC_GGO` |

## 3 Strategies

| ID | Strategy | Description |
|----|----------|-------------|
| s1_cellpose_probmap | strategy1 | Ilastik prob map -> Cellpose |
| s2_cellpose_dna_probmap | strategy2_probmap | DNA -> normalize -> gen_prob_map -> Cellpose |
| s2_cellpose_dna_direct | strategy2_direct | DNA -> normalize -> Cellpose directly |

GT comparison via `downsample_mask()` (handles 2x resolution mismatch for prob map strategies).

## Current State

- [x] Main array job 13974653 completed (91/112 CSVs)
- [ ] Retry job 13992020 in progress (21 missing tasks: 11,14,36,42-50,62,63,73,80,90,104,107-109)
- [ ] Aggregation: run `bm_aggregate.py` once all 112 CSVs exist
- [ ] Review `summary_stats.csv` and document findings

## Scripts on brb

All at `/athena/elementolab/scratch/juk4007/sc_tools/benchmark_results/v4_parallel/`:

```
bm_discover.py       # ROI manifest (done, 112 ROIs)
bm_worker.py         # per-ROI worker (SLURM_ARRAY_TASK_ID -> per_roi/task_XXXX.csv)
run_bm_array.sh      # original array job (--array=0-111)
run_bm_retry.sh      # retry (--array=11,14,36,42-50,62,63,73,80,90,104,107-109)
bm_aggregate.py      # combines CSVs -> all_results.csv + summary_stats.csv
```

## Next Steps

1. Wait for job 13992020 to finish
2. Check missing CSVs: `ssh brb 'for i in $(seq 0 111); do f=$(printf ".../task_%04d.csv" $i); [ ! -f "$f" ] && echo "$i"; done'`
3. If still missing, re-run individual workers directly (not sbatch): `ssh brb 'python bm_worker.py <id>'`
4. Run aggregation:
   ```bash
   ssh brb 'eval "$(conda shell.bash hook)" && conda activate sc_tools && python .../bm_aggregate.py'
   ```
5. Review `summary_stats.csv` — key metrics: F1, precision, recall, mean_iou per dataset per strategy

## Early Results (partial, from completed tasks)

| Strategy | F1 (approx) | Runtime/ROI |
|----------|-------------|-------------|
| s1_cellpose_probmap | 0.70-0.79 | ~450s |
| s2_cellpose_dna_probmap | ~0.81 | ~160s |
| s2_cellpose_dna_direct | TBD | TBD |

## Known Issues

- TF/numpy conflict: DeepCell/StarDist cannot run in sc_tools env (TF needs numpy<2.0, env has 2.4.2). Needs separate conda env — future work.
- Cellpose v4: Only `CellposeModel` exists; `model_type` ignored, always uses `cpsam`.
- Pip race condition: Never run `pip install -e .` inside array tasks. Already fixed in retry script (`set -e` after pip).
