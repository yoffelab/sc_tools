"""Tests for provenance infrastructure (PRV-01, PRV-02, PRV-05).

Covers:
- SHA256 streaming checksum
- Sidecar path computation, writing, reading
- Peak memory helper
- InputFile and ProvenanceRecord models
- cli_handler sidecar hook (success/error/empty-artifacts)
- h5ad uns embedding
- Leiden random_state threading and reproducibility
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


# ---------------------------------------------------------------------------
# SHA256 checksum
# ---------------------------------------------------------------------------


def test_sha256_file_correct_digest(tmp_path: Path) -> None:
    """sha256_file returns correct hex digest for known-content file."""
    import hashlib

    from sc_tools.provenance.checksum import sha256_file

    content = b"hello provenance world"
    f = tmp_path / "test.txt"
    f.write_bytes(content)

    expected = hashlib.sha256(content).hexdigest()
    assert sha256_file(f) == expected


# ---------------------------------------------------------------------------
# Sidecar path
# ---------------------------------------------------------------------------


def test_sidecar_path_for_h5ad() -> None:
    """sidecar_path_for returns .provenance.json suffix (D-02)."""
    from sc_tools.provenance.sidecar import sidecar_path_for

    result = sidecar_path_for("output.h5ad")
    assert result == Path("output.h5ad.provenance.json")


# ---------------------------------------------------------------------------
# Sidecar write / read
# ---------------------------------------------------------------------------


def test_write_sidecar_creates_json(tmp_path: Path) -> None:
    """write_sidecar creates atomic JSON file with all ProvenanceRecord fields."""
    from sc_tools.models.result import ProvenanceRecord
    from sc_tools.provenance.sidecar import read_sidecar, write_sidecar

    artifact = tmp_path / "out.h5ad"
    artifact.write_bytes(b"fake")

    record = ProvenanceRecord(
        command="preprocess run",
        params={"resolution": 0.8, "batch_key": "sample"},
        runtime_s=12.5,
        peak_memory_mb=1024.0,
    )
    sp = write_sidecar(artifact, record)
    assert sp.exists()

    data = json.loads(sp.read_text())
    assert data["command"] == "preprocess run"
    assert data["params"]["resolution"] == 0.8
    assert data["runtime_s"] == 12.5
    assert data["peak_memory_mb"] == 1024.0
    assert "sc_tools_version" in data
    assert "timestamp" in data
    assert "inputs" in data


def test_read_sidecar_returns_dict(tmp_path: Path) -> None:
    """read_sidecar returns parsed dict matching written record."""
    from sc_tools.models.result import ProvenanceRecord
    from sc_tools.provenance.sidecar import read_sidecar, write_sidecar

    artifact = tmp_path / "out.csv"
    artifact.write_text("a,b\n1,2\n")

    record = ProvenanceRecord(command="qc run", params={"modality": "visium"})
    write_sidecar(artifact, record)

    result = read_sidecar(artifact)
    assert result is not None
    assert result["command"] == "qc run"


def test_read_sidecar_returns_none_when_missing(tmp_path: Path) -> None:
    """read_sidecar returns None when no sidecar file exists."""
    from sc_tools.provenance.sidecar import read_sidecar

    assert read_sidecar(tmp_path / "nonexistent.h5ad") is None


# ---------------------------------------------------------------------------
# Peak memory
# ---------------------------------------------------------------------------


def test_get_peak_memory_mb_positive() -> None:
    """get_peak_memory_mb returns a positive float."""
    from sc_tools.provenance.sidecar import get_peak_memory_mb

    val = get_peak_memory_mb()
    assert isinstance(val, float)
    assert val > 0


# ---------------------------------------------------------------------------
# InputFile model
# ---------------------------------------------------------------------------


def test_input_file_model_fields() -> None:
    """InputFile model validates path, path_type, sha256, size_bytes (D-08)."""
    from sc_tools.models.result import InputFile

    inf = InputFile(path="data/raw.h5ad", sha256="abc123", size_bytes=1000)
    assert inf.path == "data/raw.h5ad"
    assert inf.path_type == "relative"
    assert inf.sha256 == "abc123"
    assert inf.size_bytes == 1000


# ---------------------------------------------------------------------------
# ProvenanceRecord model
# ---------------------------------------------------------------------------


def test_provenance_record_fields() -> None:
    """ProvenanceRecord contains all required fields (D-06)."""
    from sc_tools.models.result import ProvenanceRecord

    rec = ProvenanceRecord(
        command="preprocess run",
        params={"resolution": 0.8},
        runtime_s=5.0,
        peak_memory_mb=512.0,
    )
    assert rec.command == "preprocess run"
    assert rec.params == {"resolution": 0.8}
    assert isinstance(rec.inputs, list)
    assert rec.sc_tools_version is not None
    assert rec.timestamp is not None
    assert rec.runtime_s == 5.0
    assert rec.peak_memory_mb == 512.0


# ---------------------------------------------------------------------------
# cli_handler sidecar hook
# ---------------------------------------------------------------------------


def test_cli_handler_writes_sidecar_on_success(tmp_path: Path) -> None:
    """cli_handler writes sidecar when result.status=success and artifacts non-empty (D-01)."""
    from sc_tools.cli import cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status
    from sc_tools.provenance.sidecar import sidecar_path_for

    artifact = tmp_path / "output.csv"
    artifact.write_text("data")

    @cli_handler
    def fake_cmd(**kwargs):
        return CLIResult(
            status=Status.success,
            command="test cmd",
            artifacts=[str(artifact)],
            provenance=Provenance(command="test cmd"),
            message="ok",
        )

    with pytest.raises(SystemExit) as exc_info:
        fake_cmd()

    assert exc_info.value.code == 0
    sp = sidecar_path_for(artifact)
    assert sp.exists(), f"Sidecar not found at {sp}"


def test_cli_handler_no_sidecar_on_error(tmp_path: Path) -> None:
    """cli_handler does NOT write sidecar when result.status=error (D-03)."""
    from sc_tools.cli import cli_handler
    from sc_tools.errors import SCToolsUserError
    from sc_tools.provenance.sidecar import sidecar_path_for

    artifact = tmp_path / "output.csv"
    artifact.write_text("data")

    @cli_handler
    def fail_cmd(**kwargs):
        raise SCToolsUserError("bad input", suggestion="fix it")

    with pytest.raises(SystemExit):
        fail_cmd()

    sp = sidecar_path_for(artifact)
    assert not sp.exists()


def test_cli_handler_no_sidecar_when_no_artifacts() -> None:
    """cli_handler does NOT write sidecar when artifacts is empty."""
    from sc_tools.cli import cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status

    @cli_handler
    def no_artifacts_cmd(**kwargs):
        return CLIResult(
            status=Status.success,
            command="test cmd",
            artifacts=[],
            provenance=Provenance(command="test cmd"),
            message="ok",
        )

    # Should succeed without error (no sidecar to check)
    with pytest.raises(SystemExit) as exc_info:
        no_artifacts_cmd()
    assert exc_info.value.code == 0


def test_sidecar_params_include_kwargs(tmp_path: Path) -> None:
    """Sidecar params include all kwargs passed to the command (D-07)."""
    from sc_tools.cli import cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status
    from sc_tools.provenance.sidecar import read_sidecar

    artifact = tmp_path / "output.csv"
    artifact.write_text("data")

    @cli_handler
    def param_cmd(resolution=0.8, batch_key="sample", **kwargs):
        return CLIResult(
            status=Status.success,
            command="test cmd",
            artifacts=[str(artifact)],
            provenance=Provenance(command="test cmd"),
            message="ok",
        )

    with pytest.raises(SystemExit):
        param_cmd(resolution=1.5, batch_key="library_id")

    prov = read_sidecar(artifact)
    assert prov is not None
    assert prov["params"]["resolution"] == 1.5
    assert prov["params"]["batch_key"] == "library_id"


def test_input_files_convention_key(tmp_path: Path) -> None:
    """_input_files convention key in data dict is consumed by cli_handler and removed before JSON emission."""
    import io
    import sys

    from sc_tools.cli import cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status

    input_f = tmp_path / "input.h5ad"
    input_f.write_bytes(b"fake input data")
    artifact = tmp_path / "output.csv"
    artifact.write_text("data")

    @cli_handler
    def input_cmd(**kwargs):
        return CLIResult(
            status=Status.success,
            command="test cmd",
            data={"_input_files": [str(input_f)], "key": "value"},
            artifacts=[str(artifact)],
            provenance=Provenance(command="test cmd"),
            message="ok",
        )

    captured = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = captured

    try:
        with pytest.raises(SystemExit):
            input_cmd()
    finally:
        sys.stdout = old_stdout

    output = captured.getvalue()
    emitted = json.loads(output)
    # _input_files should be removed from emitted JSON
    assert "_input_files" not in emitted["data"]
    assert emitted["data"]["key"] == "value"


# ---------------------------------------------------------------------------
# h5ad uns embedding
# ---------------------------------------------------------------------------


def test_embed_provenance_in_adata(tmp_path: Path) -> None:
    """embed_provenance_in_adata writes provenance JSON to adata.uns['sct_provenance'] (D-04)."""
    h5py = pytest.importorskip("h5py")

    from sc_tools.models.result import ProvenanceRecord
    from sc_tools.provenance.sidecar import embed_provenance_in_adata

    # Create minimal h5ad-like file
    h5ad_path = tmp_path / "test.h5ad"
    with h5py.File(h5ad_path, "w") as f:
        f.create_group("uns")

    record = ProvenanceRecord(command="preprocess run", params={"resolution": 0.8})
    embed_provenance_in_adata(h5ad_path, record)

    # Verify
    with h5py.File(h5ad_path, "r") as f:
        assert "uns/sct_provenance" in f
        raw = f["uns/sct_provenance"][()]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        prov = json.loads(raw)
        assert prov["command"] == "preprocess run"


def test_write_provenance_sidecars_calls_embed_for_h5ad() -> None:
    """_write_provenance_sidecars calls embed_provenance_in_adata for h5ad artifacts (D-04)."""
    from sc_tools.cli import _write_provenance_sidecars
    from sc_tools.models.result import CLIResult, Provenance, ProvenanceRecord, Status

    result = CLIResult(
        status=Status.success,
        command="preprocess run",
        artifacts=["/tmp/output.h5ad"],
        provenance=Provenance(command="preprocess run"),
    )

    with patch("sc_tools.provenance.sidecar.embed_provenance_in_adata") as mock_embed, \
         patch("sc_tools.provenance.sidecar.write_sidecar"), \
         patch("sc_tools.provenance.sidecar.build_provenance_record") as mock_build:
        mock_build.return_value = ProvenanceRecord(command="preprocess run", params={})
        _write_provenance_sidecars(result, {}, [], 1.0, 100.0)
        mock_embed.assert_called_once()


# ---------------------------------------------------------------------------
# Leiden random_state (PRV-05)
# ---------------------------------------------------------------------------


def test_leiden_random_state_passed() -> None:
    """_leiden_cluster passes random_state to sc.tl.leiden."""
    import numpy as np

    def _fake_leiden(adata, **kwargs):
        """Mock leiden that creates the expected obs column."""
        adata.obs["leiden"] = ["0"] * adata.n_obs

    with patch("scanpy.tl.leiden", side_effect=_fake_leiden) as mock_leiden, \
         patch("scanpy.pp.neighbors"):
        from sc_tools.bm.integration import _leiden_cluster

        X = np.random.rand(20, 5)
        _leiden_cluster(X, resolution=1.0, random_state=42)
        mock_leiden.assert_called_once()
        call_kwargs = mock_leiden.call_args
        assert call_kwargs.kwargs.get("random_state") == 42 or call_kwargs[1].get("random_state") == 42


def test_leiden_random_state_default() -> None:
    """_leiden_cluster defaults to random_state=0."""
    import inspect

    from sc_tools.bm.integration import _leiden_cluster

    sig = inspect.signature(_leiden_cluster)
    assert sig.parameters["random_state"].default == 0


def test_compute_integration_metrics_has_random_state() -> None:
    """compute_integration_metrics accepts random_state parameter."""
    import inspect

    from sc_tools.bm.integration import compute_integration_metrics

    sig = inspect.signature(compute_integration_metrics)
    assert "random_state" in sig.parameters
    assert sig.parameters["random_state"].default == 0


def test_cluster_has_random_state() -> None:
    """pp.reduce.cluster accepts random_state parameter and defaults to 0."""
    import inspect

    from sc_tools.pp.reduce import cluster

    sig = inspect.signature(cluster)
    assert "random_state" in sig.parameters
    assert sig.parameters["random_state"].default == 0


def test_write_provenance_sidecars_skips_embed_for_non_h5ad() -> None:
    """_write_provenance_sidecars does NOT call embed_provenance_in_adata for non-h5ad artifacts."""
    from sc_tools.cli import _write_provenance_sidecars
    from sc_tools.models.result import CLIResult, Provenance, ProvenanceRecord, Status

    result = CLIResult(
        status=Status.success,
        command="report generate",
        artifacts=["/tmp/report.html"],
        provenance=Provenance(command="report generate"),
    )

    with patch("sc_tools.provenance.sidecar.embed_provenance_in_adata") as mock_embed, \
         patch("sc_tools.provenance.sidecar.write_sidecar"), \
         patch("sc_tools.provenance.sidecar.build_provenance_record") as mock_build:
        mock_build.return_value = ProvenanceRecord(command="report generate", params={})
        _write_provenance_sidecars(result, {}, [], 1.0, 100.0)
        mock_embed.assert_not_called()
