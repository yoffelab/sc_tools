"""Tests for sc_tools.bm strategy modules and prepare utilities."""

from __future__ import annotations

import numpy as np


class TestNormalizeIMCIntensity:
    def test_percentile(self):
        from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

        img = np.random.rand(10, 10).astype(np.float32) * 100
        result = normalize_imc_intensity(img, method="percentile")
        assert result.min() >= 0
        assert result.max() <= 1.0

    def test_arcsinh(self):
        from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

        img = np.random.rand(10, 10).astype(np.float32) * 100
        result = normalize_imc_intensity(img, method="arcsinh")
        expected = np.arcsinh(img / 5.0)
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_zscore(self):
        from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

        img = np.random.rand(10, 10).astype(np.float32) * 100
        result = normalize_imc_intensity(img, method="zscore")
        assert abs(result.mean()) < 0.5  # approximately zero-centered

    def test_uint8(self):
        from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

        img = np.random.rand(10, 10).astype(np.float32) * 100
        result = normalize_imc_intensity(img, method="uint8")
        assert result.dtype == np.uint8
        assert result.max() <= 255

    def test_multichannel(self):
        from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

        img = np.random.rand(3, 10, 10).astype(np.float32) * 100
        result = normalize_imc_intensity(img, method="percentile")
        assert result.shape == (3, 10, 10)


class TestGenerateProbabilityMap:
    def test_gaussian(self):
        from sc_tools.data.imc.benchmark.prepare import generate_probability_map

        dna = np.random.rand(50, 50).astype(np.float32) * 50
        result = generate_probability_map(dna, method="gaussian")
        assert result.shape == (50, 50, 3)
        assert result.dtype == np.float32
        # Probabilities should sum to ~1
        np.testing.assert_allclose(result.sum(axis=2), 1.0, atol=0.01)

    def test_otsu(self):
        from sc_tools.data.imc.benchmark.prepare import generate_probability_map

        dna = np.random.rand(50, 50).astype(np.float32) * 50
        result = generate_probability_map(dna, method="otsu")
        assert result.shape == (50, 50, 3)

    def test_multiscale(self):
        from sc_tools.data.imc.benchmark.prepare import generate_probability_map

        dna = np.random.rand(50, 50).astype(np.float32) * 50
        result = generate_probability_map(dna, method="multiscale")
        assert result.shape == (50, 50, 3)


class TestPostprocess:
    def test_semantic_to_instance(self):
        from sc_tools.bm.postprocess import semantic_to_instance

        # Create a binary mask with two blobs
        binary = np.zeros((50, 50), dtype=np.int32)
        binary[10:20, 10:20] = 1
        binary[30:40, 30:40] = 1

        instances = semantic_to_instance(binary)
        assert instances.dtype == np.int32
        n_cells = len(np.unique(instances)) - 1
        assert n_cells >= 2

    def test_filter_masks_by_area(self):
        from sc_tools.bm.postprocess import filter_masks_by_area

        mask = np.zeros((50, 50), dtype=np.int32)
        # Small cell (4 pixels)
        mask[0:2, 0:2] = 1
        # Large cell (100 pixels)
        mask[10:20, 10:20] = 2

        filtered = filter_masks_by_area(mask, min_area=10)
        # Small cell removed, large cell preserved (relabeled contiguously)
        assert len(np.unique(filtered)) - 1 == 1  # only 1 cell remains
        assert np.sum(filtered > 0) == 100  # large cell preserved

    def test_fill_holes(self):
        from sc_tools.bm.postprocess import fill_holes

        mask = np.zeros((50, 50), dtype=np.int32)
        mask[10:20, 10:20] = 1
        mask[14:16, 14:16] = 0  # hole

        filled = fill_holes(mask, max_hole_area=10)
        assert filled[15, 15] == 1  # hole filled


class TestDeepCellRunner:
    def test_prepare_inputs_probmap(self):
        from sc_tools.bm.deepcell_runner import _prepare_inputs

        prob_map = np.random.rand(50, 50, 3).astype(np.float32)
        nuclear, membrane = _prepare_inputs(prob_map, None, None, 1, 2)
        assert nuclear.shape == (50, 50)
        assert membrane.shape == (50, 50)

    def test_prepare_inputs_intensity(self):
        from sc_tools.bm.deepcell_runner import _prepare_inputs

        intensity = np.random.rand(5, 50, 50).astype(np.float32)
        nuclear, membrane = _prepare_inputs(intensity, [0, 1], [2, 3], 1, 2)
        assert nuclear.shape == (50, 50)
        assert membrane.shape == (50, 50)


class TestStrategyHF:
    def test_normalize_imc_to_uint8(self):
        from sc_tools.bm.strategy_hf import _normalize_imc_to_uint8

        img = np.random.rand(50, 50).astype(np.float32) * 100
        result = _normalize_imc_to_uint8(img)
        assert result.dtype == np.uint8
        assert result.max() <= 255

    def test_to_pseudo_rgb(self):
        from sc_tools.bm.strategy_hf import _to_pseudo_rgb

        gray = np.random.randint(0, 256, (50, 50), dtype=np.uint8)
        rgb = _to_pseudo_rgb(gray)
        assert rgb.shape == (50, 50, 3)
        np.testing.assert_array_equal(rgb[:, :, 0], gray)

    def test_model_registry(self):
        from sc_tools.bm.strategy_hf import HF_MODEL_REGISTRY

        assert "cellvit_256" in HF_MODEL_REGISTRY
        assert "sam_base" in HF_MODEL_REGISTRY
        assert "stardist" in HF_MODEL_REGISTRY


class TestBenchmarkConfig:
    def test_default_config(self):
        from sc_tools.data.imc.benchmark.config import BenchmarkConfig

        config = BenchmarkConfig()
        assert config.strategies == [1, 2, 3, 4]
        assert config.gpu is True
        assert config.resume is True

    def test_yaml_roundtrip(self, tmp_path):
        from sc_tools.data.imc.benchmark.config import BenchmarkConfig

        config = BenchmarkConfig(strategies=[1, 2], n_rois=50, gpu=False)
        yaml_path = tmp_path / "config.yaml"
        config.to_yaml(yaml_path)

        loaded = BenchmarkConfig.from_yaml(yaml_path)
        assert loaded.strategies == [1, 2]
        assert loaded.n_rois == 50
        assert loaded.gpu is False
