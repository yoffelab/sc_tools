"""Unit tests for sc_tools.storage."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_adata(n_obs: int = 5, n_vars: int = 4) -> ad.AnnData:
    X = np.random.rand(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame({"sample": [f"s{i}" for i in range(n_obs)]})
    var = pd.DataFrame({"gene": [f"g{i}" for i in range(n_vars)]})
    obs.index = [f"cell{i}" for i in range(n_obs)]
    var.index = [f"g{i}" for i in range(n_vars)]
    return ad.AnnData(X=X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# resolve_fs
# ---------------------------------------------------------------------------


class TestResolveFs:
    def test_local_path_str(self, tmp_path):
        from sc_tools.storage import resolve_fs

        fs, path = resolve_fs(str(tmp_path / "file.h5ad"))
        assert path.endswith("file.h5ad")

    def test_local_path_pathlike(self, tmp_path):
        from sc_tools.storage import resolve_fs

        fs, path = resolve_fs(tmp_path / "file.h5ad")
        assert "file.h5ad" in path

    def test_file_scheme(self, tmp_path):
        from sc_tools.storage import resolve_fs

        uri = (tmp_path / "x.tsv").as_uri()  # file:///...
        fs, path = resolve_fs(uri)
        assert "x.tsv" in path

    def test_missing_backend_raises_import_error(self):
        from sc_tools.storage import resolve_fs

        # Patch fsspec.url_to_fs to raise ImportError (simulates missing backend)
        with patch("fsspec.url_to_fs", side_effect=ImportError("no s3fs")):
            with pytest.raises(ImportError, match="Storage backend"):
                resolve_fs("s3://bucket/key")


# ---------------------------------------------------------------------------
# open_file
# ---------------------------------------------------------------------------


class TestOpenFile:
    def test_read_local_text(self, tmp_path):
        from sc_tools.storage import open_file

        f = tmp_path / "data.txt"
        f.write_text("hello\nworld\n")
        with open_file(str(f), "r") as fh:
            content = fh.read()
        assert "hello" in content

    def test_write_local_bytes(self, tmp_path):
        from sc_tools.storage import open_file

        f = tmp_path / "out.bin"
        with open_file(str(f), "wb") as fh:
            fh.write(b"\x00\x01\x02")
        assert f.read_bytes() == b"\x00\x01\x02"


# ---------------------------------------------------------------------------
# with_local_copy
# ---------------------------------------------------------------------------


class TestWithLocalCopy:
    def test_local_path_noop(self, tmp_path):
        from sc_tools.storage import with_local_copy

        f = tmp_path / "file.h5ad"
        f.write_bytes(b"dummy")
        with with_local_copy(str(f)) as local:
            assert local.exists()
            assert local == f or str(local) == str(f)

    def test_remote_downloads_to_tmp(self, tmp_path):
        from sc_tools.storage import with_local_copy

        # Patch fsspec so "remote" URI resolves to a mock FS
        content = b"fake_h5ad_bytes"
        mock_fs = MagicMock()
        mock_fs.get = lambda src, dst: Path(dst).write_bytes(content)

        def fake_url_to_fs(uri):
            return mock_fs, "/remote/path/file.h5ad"

        with (
            patch("fsspec.url_to_fs", side_effect=fake_url_to_fs),
            patch(
                "sc_tools.storage._is_local",
                side_effect=lambda fs: False,
            ),
        ):
            with with_local_copy("sftp://host//remote/path/file.h5ad") as local:
                assert local.suffix == ".h5ad"
                assert local.read_bytes() == content
            # tmp file is cleaned up after context
            assert not local.exists()


# ---------------------------------------------------------------------------
# smart_read_h5ad
# ---------------------------------------------------------------------------


class TestSmartReadH5ad:
    def test_local_h5ad(self, tmp_path):
        from sc_tools.storage import smart_read_h5ad

        adata = _make_adata()
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        loaded = smart_read_h5ad(str(path))
        assert loaded.n_obs == adata.n_obs
        assert loaded.n_vars == adata.n_vars

    def test_local_pathlike(self, tmp_path):
        from sc_tools.storage import smart_read_h5ad

        adata = _make_adata()
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        loaded = smart_read_h5ad(path)
        assert loaded.n_obs == adata.n_obs


# ---------------------------------------------------------------------------
# smart_write_checkpoint
# ---------------------------------------------------------------------------


class TestSmartWriteCheckpoint:
    def test_local_h5ad(self, tmp_path):
        from sc_tools.storage import smart_read_h5ad, smart_write_checkpoint

        adata = _make_adata()
        uri = str(tmp_path / "out.h5ad")
        smart_write_checkpoint(adata, uri, fmt="h5ad")

        assert Path(uri).exists()
        loaded = smart_read_h5ad(uri)
        assert loaded.n_obs == adata.n_obs

    def test_creates_parent_dir(self, tmp_path):
        from sc_tools.storage import smart_write_checkpoint

        adata = _make_adata()
        uri = str(tmp_path / "deep" / "nested" / "out.h5ad")
        smart_write_checkpoint(adata, uri)
        assert Path(uri).exists()

    def test_unknown_format_raises(self, tmp_path):
        from sc_tools.storage import smart_write_checkpoint

        adata = _make_adata()
        with pytest.raises(ValueError, match="Unknown checkpoint format"):
            smart_write_checkpoint(adata, str(tmp_path / "out.xyz"), fmt="xyz")

    def test_zarr_local(self, tmp_path):
        pytest.importorskip("zarr")
        from sc_tools.storage import smart_write_checkpoint

        adata = _make_adata()
        uri = str(tmp_path / "out.zarr")
        smart_write_checkpoint(adata, uri, fmt="zarr")
        assert Path(uri).exists()


# ---------------------------------------------------------------------------
# smart_read_csv
# ---------------------------------------------------------------------------


class TestSmartReadCsv:
    def test_local_csv(self, tmp_path):
        from sc_tools.storage import smart_read_csv

        f = tmp_path / "data.csv"
        f.write_text("a,b,c\n1,2,3\n4,5,6\n")
        df = smart_read_csv(str(f))
        assert list(df.columns) == ["a", "b", "c"]
        assert len(df) == 2

    def test_local_tsv(self, tmp_path):
        from sc_tools.storage import smart_read_csv

        f = tmp_path / "data.tsv"
        f.write_text("a\tb\n1\t2\n3\t4\n")
        df = smart_read_csv(str(f), sep="\t")
        assert list(df.columns) == ["a", "b"]
        assert len(df) == 2
