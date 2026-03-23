"""Tests for provenance CLI commands and lineage trace engine.

Tests cover:
- trace_lineage: BFS walk with cycle detection, adata.uns fallback, SHA256 relocation
- CLI commands: sct provenance show, sct provenance trace
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# trace_lineage tests
# ---------------------------------------------------------------------------


class TestTraceLineageNoSidecar:
    """trace_lineage on file with no sidecar returns origin step."""

    def test_no_sidecar_returns_single_origin_step(self, tmp_path: Path) -> None:
        from sc_tools.provenance.trace import trace_lineage

        f = tmp_path / "raw_data.csv"
        f.write_text("data")

        steps = trace_lineage(str(f), project_dir=str(tmp_path))

        assert len(steps) == 1
        assert steps[0]["file"] == str(f)
        assert steps[0]["command"] is None
        assert steps[0]["timestamp"] is None
        assert steps[0]["note"] == "origin (no provenance)"


class TestTraceLineageWithSidecar:
    """trace_lineage on file with sidecar returns step with command info."""

    def test_single_sidecar_returns_command_info(self, tmp_path: Path) -> None:
        from sc_tools.provenance.trace import trace_lineage

        f = tmp_path / "output.h5ad"
        f.write_bytes(b"fake h5ad")
        sidecar = Path(str(f) + ".provenance.json")
        sidecar.write_text(
            json.dumps(
                {
                    "command": "preprocess run",
                    "timestamp": "2026-03-23T10:00:00Z",
                    "params": {},
                    "inputs": [],
                }
            )
        )

        steps = trace_lineage(str(f), project_dir=str(tmp_path))

        assert len(steps) == 1
        assert steps[0]["command"] == "preprocess run"
        assert steps[0]["timestamp"] == "2026-03-23T10:00:00Z"
        assert steps[0]["note"] is None


class TestTraceLineageChain:
    """trace_lineage follows input references recursively."""

    def test_three_step_chain_produces_chronological_output(self, tmp_path: Path) -> None:
        """A -> B -> C chain: C references B, B references A. Output oldest first."""
        from sc_tools.provenance.trace import trace_lineage

        # Create files
        file_a = tmp_path / "raw.csv"
        file_a.write_text("raw data")

        file_b = tmp_path / "filtered.h5ad"
        file_b.write_bytes(b"filtered")

        file_c = tmp_path / "normalized.h5ad"
        file_c.write_bytes(b"normalized")

        # Sidecar for B: references A as input
        sidecar_b = Path(str(file_b) + ".provenance.json")
        sidecar_b.write_text(
            json.dumps(
                {
                    "command": "qc run",
                    "timestamp": "2026-03-23T09:00:00Z",
                    "params": {},
                    "inputs": [{"path": str(file_a), "sha256": "abc", "size_bytes": 8, "path_type": "absolute"}],
                }
            )
        )

        # Sidecar for C: references B as input
        sidecar_c = Path(str(file_c) + ".provenance.json")
        sidecar_c.write_text(
            json.dumps(
                {
                    "command": "preprocess run",
                    "timestamp": "2026-03-23T10:00:00Z",
                    "params": {},
                    "inputs": [{"path": str(file_b), "sha256": "def", "size_bytes": 8, "path_type": "absolute"}],
                }
            )
        )

        steps = trace_lineage(str(file_c), project_dir=str(tmp_path))

        # 3 steps: A (origin), B (qc run), C (preprocess run)
        assert len(steps) == 3
        # Chronological order: oldest first (A, then B, then C)
        assert steps[0]["file"] == str(file_a)
        assert steps[0]["note"] == "origin (no provenance)"
        assert steps[1]["file"] == str(file_b)
        assert steps[1]["command"] == "qc run"
        assert steps[2]["file"] == str(file_c)
        assert steps[2]["command"] == "preprocess run"


class TestTraceLineageCycleDetection:
    """trace_lineage handles cycles without infinite loop."""

    def test_cycle_detection_completes(self, tmp_path: Path) -> None:
        """A references B, B references A -- should not hang."""
        from sc_tools.provenance.trace import trace_lineage

        file_a = tmp_path / "a.h5ad"
        file_a.write_bytes(b"a data")

        file_b = tmp_path / "b.h5ad"
        file_b.write_bytes(b"b data")

        # A references B
        sidecar_a = Path(str(file_a) + ".provenance.json")
        sidecar_a.write_text(
            json.dumps(
                {
                    "command": "step1",
                    "timestamp": "2026-03-23T09:00:00Z",
                    "params": {},
                    "inputs": [{"path": str(file_b), "sha256": "x", "size_bytes": 6, "path_type": "absolute"}],
                }
            )
        )

        # B references A
        sidecar_b = Path(str(file_b) + ".provenance.json")
        sidecar_b.write_text(
            json.dumps(
                {
                    "command": "step2",
                    "timestamp": "2026-03-23T08:00:00Z",
                    "params": {},
                    "inputs": [{"path": str(file_a), "sha256": "y", "size_bytes": 6, "path_type": "absolute"}],
                }
            )
        )

        steps = trace_lineage(str(file_a), project_dir=str(tmp_path))
        assert len(steps) == 2


class TestTraceLineageUnsfallback:
    """trace_lineage falls back to adata.uns when sidecar missing for h5ad."""

    def test_reads_provenance_from_uns(self, tmp_path: Path) -> None:
        h5py = pytest.importorskip("h5py")

        h5ad_path = tmp_path / "embedded.h5ad"
        prov_data = {
            "command": "qc run",
            "timestamp": "2026-03-23T09:00:00Z",
            "params": {"min_genes": 200},
            "inputs": [],
        }

        # Write minimal h5ad with uns/sct_provenance
        with h5py.File(h5ad_path, "w") as f:
            uns = f.create_group("uns")
            uns.create_dataset("sct_provenance", data=json.dumps(prov_data))

        from sc_tools.provenance.trace import trace_lineage

        steps = trace_lineage(str(h5ad_path), project_dir=str(tmp_path))

        assert len(steps) == 1
        assert steps[0]["command"] == "qc run"
        assert steps[0]["timestamp"] == "2026-03-23T09:00:00Z"


class TestTraceLineageMissingIntermediate:
    """trace_lineage handles missing intermediate files gracefully."""

    def test_missing_input_marked_as_origin(self, tmp_path: Path) -> None:
        from sc_tools.provenance.trace import trace_lineage

        output = tmp_path / "result.h5ad"
        output.write_bytes(b"result")

        # Sidecar references a file that does not exist
        sidecar = Path(str(output) + ".provenance.json")
        sidecar.write_text(
            json.dumps(
                {
                    "command": "preprocess run",
                    "timestamp": "2026-03-23T10:00:00Z",
                    "params": {},
                    "inputs": [
                        {
                            "path": str(tmp_path / "nonexistent.h5ad"),
                            "sha256": "deadbeef",
                            "size_bytes": 100,
                            "path_type": "absolute",
                        }
                    ],
                }
            )
        )

        steps = trace_lineage(str(output), project_dir=str(tmp_path))

        assert len(steps) == 2
        # The missing file should be marked as origin
        assert steps[0]["note"] == "origin (no provenance)"
        assert steps[1]["command"] == "preprocess run"


class TestReadUnsProvenance:
    """_read_uns_provenance reads provenance from h5ad without full AnnData load."""

    def test_reads_uns_json(self, tmp_path: Path) -> None:
        h5py = pytest.importorskip("h5py")

        h5ad_path = tmp_path / "test.h5ad"
        prov = {"command": "normalize", "timestamp": "2026-01-01T00:00:00Z", "params": {}, "inputs": []}

        with h5py.File(h5ad_path, "w") as f:
            uns = f.create_group("uns")
            uns.create_dataset("sct_provenance", data=json.dumps(prov))

        from sc_tools.provenance.trace import _read_uns_provenance

        result = _read_uns_provenance(Path(h5ad_path))
        assert result is not None
        assert result["command"] == "normalize"

    def test_returns_none_when_no_key(self, tmp_path: Path) -> None:
        h5py = pytest.importorskip("h5py")

        h5ad_path = tmp_path / "no_prov.h5ad"
        with h5py.File(h5ad_path, "w") as f:
            f.create_group("uns")

        from sc_tools.provenance.trace import _read_uns_provenance

        result = _read_uns_provenance(Path(h5ad_path))
        assert result is None


class TestSHA256Relocation:
    """SHA256-based relocation finds file at new path when original is missing."""

    def test_finds_relocated_file_by_sha256(self, tmp_path: Path) -> None:
        from sc_tools.provenance.trace import trace_lineage

        # Create the output file
        output = tmp_path / "result.csv"
        output.write_text("result data")

        # Create the input file at a new location (not the path in the sidecar)
        new_dir = tmp_path / "moved"
        new_dir.mkdir()
        relocated = new_dir / "input_moved.csv"
        content = b"original input data"
        relocated.write_bytes(content)
        sha = hashlib.sha256(content).hexdigest()

        # The sidecar references the old (missing) path but has the correct SHA256
        old_path = str(tmp_path / "original_input.csv")
        sidecar = Path(str(output) + ".provenance.json")
        sidecar.write_text(
            json.dumps(
                {
                    "command": "process",
                    "timestamp": "2026-03-23T10:00:00Z",
                    "params": {},
                    "inputs": [
                        {
                            "path": old_path,
                            "sha256": sha,
                            "size_bytes": len(content),
                            "path_type": "absolute",
                        }
                    ],
                }
            )
        )

        steps = trace_lineage(str(output), project_dir=str(tmp_path))

        assert len(steps) == 2
        # The relocated file should be found and marked with relocation note
        relocated_step = steps[0]
        assert "relocated" in relocated_step["note"]
        # The main output step
        assert steps[1]["command"] == "process"
