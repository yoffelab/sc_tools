"""Tests for sc_tools.bm.runner — benchmark orchestration."""

from __future__ import annotations

import pandas as pd


class TestAggregateResults:
    def test_aggregate_empty(self):
        from sc_tools.bm.runner import aggregate_results

        result = aggregate_results(pd.DataFrame())
        assert all(len(v) == 0 for v in result.values())

    def test_aggregate_groups(self):
        from sc_tools.bm.runner import aggregate_results

        df = pd.DataFrame(
            {
                "strategy": [1, 1, 2, 2],
                "method": ["cellpose", "cellpose", "stardist", "stardist"],
                "dataset": ["ds1", "ds2", "ds1", "ds2"],
                "tissue": ["lung", "lung", "lung", "colon"],
                "n_cells": [100, 150, 80, 120],
                "boundary_regularity": [0.8, 0.85, 0.7, 0.75],
            }
        )

        result = aggregate_results(df)
        assert "by_method" in result
        assert "by_dataset" in result
        assert "by_tissue" in result
        assert len(result["by_method"]) == 2
        assert len(result["by_dataset"]) == 4


class TestSLURM:
    def test_generate_sbatch(self, tmp_path):
        from sc_tools.bm.slurm import generate_benchmark_sbatch
        from sc_tools.data.imc.benchmark.config import BenchmarkConfig

        catalog = pd.DataFrame(
            {
                "dataset": ["ds1"] * 5,
                "roi_id": [f"roi{i}" for i in range(5)],
                "tiff_path": [f"/tmp/roi{i}.tiff" for i in range(5)],
            }
        )
        config = BenchmarkConfig(strategies=[1, 2], slurm_chunk_size=3)

        scripts = generate_benchmark_sbatch(
            catalog, config, output_dir=tmp_path, conda_env="test_env"
        )
        assert len(scripts) > 0
        for s in scripts:
            assert s.is_file()
            content = s.read_text()
            assert "sbatch" in content.lower() or "SBATCH" in content

    def test_collect_no_results(self, tmp_path):
        from sc_tools.bm.slurm import collect_benchmark_results

        df = collect_benchmark_results(tmp_path)
        assert len(df) == 0
