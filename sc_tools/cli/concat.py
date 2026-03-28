"""Sample concatenation command (CONCAT-01, CONCAT-02, CONCAT-03).

Provides ``sct concat`` for merging same-modality h5ad files with
spatial metadata preservation and provenance tracking.
"""

from __future__ import annotations

import logging

import typer

logger = logging.getLogger(__name__)


def register_concat(app: typer.Typer) -> None:
    """Register the concat command on the given Typer app."""
    from sc_tools.cli import cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status

    @app.command("concat")
    @cli_handler
    def concat(
        input: list[str] = typer.Option(  # noqa: A002
            ...,
            "--input",
            "-i",
            help="Input h5ad file paths (at least 2 required)",
        ),
        output: str = typer.Option(..., "--output", "-o", help="Output merged h5ad path"),
        batch_key: str = typer.Option("sample", "--batch-key", help="Obs column for sample identity"),
        dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs without merging"),
        force: bool = typer.Option(False, "--force", help="Bypass memory safety guard"),
    ) -> None:
        """Concatenate same-modality h5ad files into a single AnnData."""
        from pathlib import Path

        import anndata as ad
        import numpy as np

        from sc_tools.errors import SCToolsDataError, SCToolsUserError
        from sc_tools.storage import smart_write_checkpoint

        input_paths = [Path(p) for p in input]

        # Validate at least 2 inputs
        if len(input_paths) < 2:
            raise SCToolsUserError(
                f"At least 2 input files required, got {len(input_paths)}",
                suggestion="Provide 2 or more --input/-i paths",
            )

        # Validate all input paths exist
        missing = [str(p) for p in input_paths if not p.exists()]
        if missing:
            raise SCToolsUserError(
                f"Input files not found: {', '.join(missing)}",
                suggestion="Check file paths and ensure files exist",
            )

        # Dry run: validate inputs only, no merge or write
        if dry_run:
            return CLIResult(
                status=Status.success,
                command="concat",
                data={
                    "dry_run": True,
                    "n_inputs": len(input_paths),
                    "inputs": [str(p) for p in input_paths],
                },
                provenance=Provenance(command="concat"),
                message=f"Dry run: {len(input_paths)} inputs validated, no output written",
            )

        # Pre-flight: check gene overlap via h5py (lightweight)
        try:
            import h5py

            all_var_names = []
            for p in input_paths:
                with h5py.File(p, "r") as f:
                    if "var" in f and "_index" in f["var"]:
                        vn = [v.decode() if isinstance(v, bytes) else v for v in f["var"]["_index"][:]]
                        all_var_names.append(set(vn))

            if len(all_var_names) >= 2:
                intersection = all_var_names[0]
                union = all_var_names[0]
                for vn in all_var_names[1:]:
                    intersection = intersection & vn
                    union = union | vn
                if union:
                    overlap = len(intersection) / len(union)
                    if overlap < 0.8:
                        logger.warning(
                            "Gene overlap across inputs is %.1f%% (< 80%%). "
                            "Consider verifying inputs are same modality.",
                            overlap * 100,
                        )
        except Exception:
            logger.debug("Pre-flight gene overlap check skipped", exc_info=True)

        # Load all inputs
        adatas = []
        sample_names = []
        for p in input_paths:
            adata = ad.read_h5ad(p)
            # Cast categorical obs and var columns to str to prevent
            # duplicate-category errors during ad.concat across samples
            # with differing category sets.
            for col in adata.obs.columns:
                if hasattr(adata.obs[col], "cat"):
                    adata.obs[col] = adata.obs[col].astype(str)
            for col in adata.var.columns:
                if hasattr(adata.var[col], "cat"):
                    adata.var[col] = adata.var[col].astype(str)
            adatas.append(adata)
            # Use parent directory name as sample key (more unique than file stem,
            # which is often "adata.008um" for all Visium HD inputs).
            sample_names.append(p.parent.name)

        # Collect spatial metadata before concat — uns_merge="unique" drops
        # entries that differ across samples (i.e. every sample's images).
        # We rebuild uns["spatial"] manually after concat.
        spatial_metadata: dict = {}
        for adata_i in adatas:
            for lib_id, lib_data in adata_i.uns.get("spatial", {}).items():
                spatial_metadata[lib_id] = lib_data

        # Concatenate (uns_merge="unique" handles non-spatial uns keys)
        merged = ad.concat(
            adatas,
            join="outer",
            uns_merge="unique",
            index_unique="-",
            label=batch_key,
            keys=sample_names,
        )

        # Restore full spatial metadata with all library IDs and images
        if spatial_metadata:
            merged.uns["spatial"] = spatial_metadata

        # Ensure spatial coordinates are numpy array
        if "spatial" in merged.obsm:
            merged.obsm["spatial"] = np.array(merged.obsm["spatial"])

        # Post-concat verification: check all spatial keys present
        merged_spatial = merged.uns.get("spatial", {})
        for adata_i, name in zip(adatas, sample_names, strict=True):
            expected_keys = set(adata_i.uns.get("spatial", {}).keys())
            missing_keys = expected_keys - set(merged_spatial.keys())
            if missing_keys:
                raise SCToolsDataError(
                    f"Spatial keys lost during concat for {name}: {missing_keys}",
                    suggestion="This may indicate an anndata version incompatibility",
                )

        # Convert any Arrow- or extension-backed arrays to numpy-compatible dtypes
        # before writing — anndata's h5py backend cannot serialize ArrowStringArray
        # or other pandas extension types (pd.StringDtype with pyarrow storage).
        import numpy as np
        import pandas as pd

        # Disable pyarrow string inference so pd.Index / pd.Categorical / column
        # assignments don't silently re-infer ArrowStringArray.
        _infer_string_orig = pd.options.future.infer_string
        try:
            pd.options.future.infer_string = False

            # Convert all string columns to Categoricals first so anndata's
            # internal strings_to_categoricals() call during write_h5ad has
            # nothing left to convert (which would produce Arrow-backed categories).
            merged.strings_to_categoricals()

            def _coerce_df_to_numpy(df: pd.DataFrame) -> None:
                # Coerce all Categorical columns to use numpy object arrays for
                # both values and categories — anndata cannot write Arrow-backed
                # categories to h5py.
                for col in df.columns:
                    dtype = df[col].dtype
                    if isinstance(dtype, pd.CategoricalDtype):
                        new_cats = pd.Index(dtype.categories.tolist(), dtype=object)
                        new_vals = np.array(df[col].astype(object).tolist(), dtype=object)
                        df[col] = pd.Categorical(new_vals, categories=new_cats)
                    elif not isinstance(dtype, np.dtype):
                        df[col] = pd.array(df[col].tolist(), dtype=object)
                if not isinstance(df.index.dtype, np.dtype):
                    df.index = pd.Index(df.index.tolist(), dtype=object)

            _coerce_df_to_numpy(merged.obs)
            _coerce_df_to_numpy(merged.var)
        finally:
            pd.options.future.infer_string = _infer_string_orig

        # Write output
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        smart_write_checkpoint(merged, output_path)

        return CLIResult(
            status=Status.success,
            command="concat",
            data={
                "n_samples": len(adatas),
                "n_obs": merged.n_obs,
                "n_vars": merged.n_vars,
                "spatial_keys_preserved": [str(k) for k in merged.uns.get("spatial", {}).keys()],
                "_input_files": [str(p) for p in input_paths],
            },
            artifacts=[str(output_path)],
            provenance=Provenance(command="concat"),
            message=f"Concatenated {len(adatas)} samples -> {merged.n_obs} obs x {merged.n_vars} vars",
        )
