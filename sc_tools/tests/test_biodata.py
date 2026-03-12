"""Unit tests for sc_tools.biodata platform registry."""

from __future__ import annotations

import pytest

from sc_tools.biodata import (
    KNOWN_PLATFORMS,
    PlatformSpec,
    get_modality_for_platform,
    get_platform,
    list_modalities,
    list_platforms,
    list_platforms_by_modality,
    platform_for_project,
    register_platform,
)


class TestPlatformSpec:
    def test_dataclass_fields(self):
        spec = PlatformSpec(
            name="test",
            label="Test Platform",
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="droplet",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor="test_vendor",
        )
        assert spec.name == "test"
        assert spec.spatial is False
        assert spec.defaults == {}

    def test_frozen(self):
        spec = get_platform("xenium")
        with pytest.raises(AttributeError):
            spec.name = "changed"  # type: ignore[misc]


class TestKnownPlatforms:
    def test_minimum_count(self):
        """Plan specifies 60+ platforms; verify we have at least that."""
        assert len(KNOWN_PLATFORMS) >= 60

    def test_all_categories_present(self):
        categories = {s.category for s in KNOWN_PLATFORMS.values()}
        assert "spatial_seq" in categories
        assert "image" in categories
        assert "rnaseq" in categories
        assert "epigenomics" in categories
        assert "genome_seq" in categories

    def test_key_platforms_exist(self):
        expected = [
            "visium",
            "visium_hd",
            "visium_hd_cell",
            "xenium",
            "cosmx_1k",
            "cosmx_6k",
            "imc",
            "mibi",
            "chromium_3p",
            "chromium_5p",
            "merscope",
            "atac_seq",
            "illumina_wgs",
            "pacbio_hifi",
            "phenocycler",
            "he",
            "cite_seq",
        ]
        for name in expected:
            assert name in KNOWN_PLATFORMS, f"Missing platform: {name}"

    def test_spatial_flag_consistency(self):
        for name, spec in KNOWN_PLATFORMS.items():
            if spec.category == "spatial_seq":
                assert spec.spatial is True, f"{name} is spatial_seq but spatial=False"

    def test_all_specs_have_required_fields(self):
        for name, spec in KNOWN_PLATFORMS.items():
            assert spec.name == name, f"Key '{name}' != spec.name '{spec.name}'"
            assert spec.label, f"{name} missing label"
            assert spec.biodata_type, f"{name} missing biodata_type"
            assert spec.category, f"{name} missing category"
            assert spec.measurement, f"{name} missing measurement"
            assert spec.resolution, f"{name} missing resolution"
            assert spec.vendor, f"{name} missing vendor"


class TestGetPlatform:
    def test_known_platform(self):
        spec = get_platform("visium")
        assert spec.name == "visium"
        assert spec.biodata_type == "spatial_seq"
        assert spec.vendor == "10x_genomics"

    def test_unknown_platform_raises(self):
        with pytest.raises(KeyError, match="Unknown platform"):
            get_platform("nonexistent_platform")


class TestListPlatforms:
    def test_no_filters(self):
        result = list_platforms()
        assert len(result) == len(KNOWN_PLATFORMS)

    def test_filter_by_category(self):
        spatial = list_platforms(category="spatial_seq")
        assert all(s.category == "spatial_seq" for s in spatial)
        assert len(spatial) > 10  # plenty of spatial platforms

    def test_filter_by_vendor(self):
        tenx = list_platforms(vendor="10x_genomics")
        assert all(s.vendor == "10x_genomics" for s in tenx)
        assert len(tenx) >= 5

    def test_filter_by_spatial(self):
        spatial = list_platforms(spatial=True)
        assert all(s.spatial is True for s in spatial)
        non_spatial = list_platforms(spatial=False)
        assert all(s.spatial is False for s in non_spatial)
        assert len(spatial) + len(non_spatial) == len(KNOWN_PLATFORMS)

    def test_filter_by_measurement(self):
        protein = list_platforms(measurement="protein")
        assert all(s.measurement == "protein" for s in protein)
        assert len(protein) >= 5

    def test_filter_by_resolution(self):
        bulk = list_platforms(resolution="bulk")
        assert all(s.resolution == "bulk" for s in bulk)

    def test_combined_filters(self):
        result = list_platforms(category="image", measurement="protein")
        assert all(s.category == "image" and s.measurement == "protein" for s in result)

    def test_sorted_by_name(self):
        result = list_platforms()
        names = [s.name for s in result]
        assert names == sorted(names)


class TestRegisterPlatform:
    def test_register_new(self):
        custom = PlatformSpec(
            name="custom_test_platform",
            label="Custom Test",
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="custom",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor="test",
        )
        register_platform("custom_test_platform", custom)
        assert get_platform("custom_test_platform") == custom
        # Cleanup
        del KNOWN_PLATFORMS["custom_test_platform"]

    def test_overwrite_existing(self):
        original = get_platform("visium")
        custom = PlatformSpec(
            name="visium",
            label="Custom Visium",
            biodata_type="spatial_seq",
            category="spatial_seq",
            subcategory="custom",
            measurement="rna",
            resolution="spot",
            spatial=True,
            vendor="custom",
        )
        register_platform("visium", custom)
        assert get_platform("visium").label == "Custom Visium"
        # Restore
        register_platform("visium", original)


class TestPlatformForProject:
    def test_known_project_platform(self):
        spec = platform_for_project("visium")
        assert spec is not None
        assert spec.name == "visium"

    def test_imc_project_platform(self):
        spec = platform_for_project("imc")
        assert spec is not None
        assert spec.measurement == "protein"

    def test_unknown_returns_none(self):
        assert platform_for_project("totally_unknown") is None


class TestModalityField:
    def test_all_platforms_have_modality(self):
        """Every registered platform should have a non-empty modality string."""
        for name, spec in KNOWN_PLATFORMS.items():
            assert spec.modality, f"Platform '{name}' has empty modality"

    def test_imc_modality(self):
        spec = get_platform("imc")
        assert spec.modality == "Spatial Proteomics - Mass Spec"

    def test_visium_modality(self):
        spec = get_platform("visium")
        assert spec.modality == "Spatial Transcriptomics - Sequencing"

    def test_xenium_modality(self):
        spec = get_platform("xenium")
        assert spec.modality == "Spatial Transcriptomics - ISH"

    def test_chromium_modality(self):
        spec = get_platform("chromium_3p")
        assert spec.modality == "Single-Cell RNA-seq - Droplet"

    def test_illumina_wgs_modality(self):
        spec = get_platform("illumina_wgs")
        assert spec.modality == "Genome Sequencing - Short Read"

    def test_atac_seq_modality(self):
        spec = get_platform("atac_seq")
        assert spec.modality == "Epigenomics - Bulk"

    def test_he_modality(self):
        spec = get_platform("he")
        assert spec.modality == "Histology"

    def test_cite_seq_modality(self):
        spec = get_platform("cite_seq")
        assert spec.modality == "Multiomics"

    def test_geomx_modality(self):
        spec = get_platform("geomx")
        assert spec.modality == "Spatial Transcriptomics - Region Capture"

    def test_phenocycler_modality(self):
        spec = get_platform("phenocycler")
        assert spec.modality == "Spatial Proteomics - Cyclic IF"


class TestPlatformVersion:
    def test_default_empty_string(self):
        spec = get_platform("visium")
        assert spec.platform_version == ""

    def test_custom_platform_with_version(self):
        custom = PlatformSpec(
            name="test_v2",
            label="Test v2",
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="custom",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor="test",
            modality="Test Modality",
            platform_version="v2",
        )
        assert custom.platform_version == "v2"


class TestPlatformDefaults:
    def test_imc_defaults(self):
        spec = get_platform("imc")
        assert spec.defaults.get("staining_protocol") == "IMC"
        assert spec.defaults.get("image_type") == "multiplexed"

    def test_mibi_defaults(self):
        spec = get_platform("mibi")
        assert spec.defaults.get("staining_protocol") == "MIBI"
        assert spec.defaults.get("image_type") == "multiplexed"

    def test_phenocycler_defaults(self):
        spec = get_platform("phenocycler")
        assert spec.defaults.get("staining_protocol") == "CODEX"
        assert spec.defaults.get("image_type") == "multiplexed"

    def test_he_defaults(self):
        spec = get_platform("he")
        assert spec.defaults.get("staining_protocol") == "H&E"
        assert spec.defaults.get("image_type") == "he"

    def test_visium_hd_defaults(self):
        spec = get_platform("visium_hd")
        assert spec.defaults.get("bin_size_um") == 8.0
        assert spec.defaults.get("spatial_resolution") == "spot"

    def test_xenium_defaults(self):
        spec = get_platform("xenium")
        assert spec.defaults.get("spatial_resolution") == "single_cell"
        assert spec.defaults.get("coordinate_system") == "micron"

    def test_chromium_3p_defaults(self):
        spec = get_platform("chromium_3p")
        assert spec.defaults.get("chemistry") == "chromium_v3"
        assert spec.defaults.get("library_type") == "single_cell"

    def test_illumina_wgs_defaults(self):
        spec = get_platform("illumina_wgs")
        assert spec.defaults.get("sequencing_type") == "wgs"
        assert "illumina" in spec.defaults.get("sequencing_platform", "")

    def test_pacbio_hifi_defaults(self):
        spec = get_platform("pacbio_hifi")
        assert spec.defaults.get("sequencing_type") == "wgs"
        assert "pacbio" in spec.defaults.get("sequencing_platform", "")


class TestListModalities:
    def test_returns_sorted_list(self):
        result = list_modalities()
        assert result == sorted(result)
        assert len(result) >= 15

    def test_filter_by_biodata_type(self):
        spatial_seq_mods = list_modalities(biodata_type="spatial_seq")
        assert "Spatial Transcriptomics - Sequencing" in spatial_seq_mods
        assert "Spatial Transcriptomics - ISH" in spatial_seq_mods
        assert "Spatial Transcriptomics - Region Capture" in spatial_seq_mods
        # Should NOT include non-spatial_seq modalities
        assert "Histology" not in spatial_seq_mods
        assert "Bulk RNA-seq" not in spatial_seq_mods

    def test_filter_by_image(self):
        image_mods = list_modalities(biodata_type="image")
        assert "Spatial Proteomics - Mass Spec" in image_mods
        assert "Histology" in image_mods

    def test_filter_by_rnaseq(self):
        rna_mods = list_modalities(biodata_type="rnaseq")
        assert "Single-Cell RNA-seq - Droplet" in rna_mods
        assert "Bulk RNA-seq" in rna_mods
        assert "Multiomics" in rna_mods

    def test_no_empty_strings(self):
        result = list_modalities()
        assert "" not in result


class TestListPlatformsByModality:
    def test_mass_spec_imaging(self):
        result = list_platforms_by_modality("Spatial Proteomics - Mass Spec")
        names = [s.name for s in result]
        assert "imc" in names
        assert "mibi" in names
        assert "maldi_ims" in names
        assert len(result) == 3

    def test_droplet_scrna(self):
        result = list_platforms_by_modality("Single-Cell RNA-seq - Droplet")
        names = [s.name for s in result]
        assert "chromium_3p" in names
        assert "dropseq" in names

    def test_empty_for_unknown(self):
        result = list_platforms_by_modality("Nonexistent Modality")
        assert result == []

    def test_sorted_by_name(self):
        result = list_platforms_by_modality("Spatial Transcriptomics - ISH")
        names = [s.name for s in result]
        assert names == sorted(names)


class TestGetModalityForPlatform:
    def test_known(self):
        assert get_modality_for_platform("imc") == "Spatial Proteomics - Mass Spec"
        assert get_modality_for_platform("visium") == "Spatial Transcriptomics - Sequencing"

    def test_unknown_raises(self):
        with pytest.raises(KeyError):
            get_modality_for_platform("nonexistent_platform")
