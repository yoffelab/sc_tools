"""MCP server: sc-tools analysis tools.

Exposes sc_tools capabilities as callable MCP tools so Claude Code can
invoke analysis steps directly.

Tools
-----
    validate_checkpoint         -- validate an AnnData checkpoint
    generate_qc_report          -- generate a QC HTML report
    score_gene_signatures       -- score gene signatures on an AnnData
    run_deconvolution           -- run cell-type deconvolution
    generate_sbatch_spaceranger -- generate a Space Ranger SLURM script
    generate_sbatch_imc         -- generate an IMC pipeline SLURM script
    collect_batch_manifests     -- collect per-batch TSVs into all_samples.tsv
    load_sample                 -- load one sample into AnnData (local path)

Start the server::

    python -m sc_tools.mcp.tools_server

Requires: pip install "sc-tools[mcp]"
"""

from __future__ import annotations

import logging

from sc_tools.models.result import CLIResult, Provenance, Status

logger = logging.getLogger(__name__)

try:
    from mcp.server.fastmcp import FastMCP
except ImportError as _mcp_err:
    raise ImportError(
        "mcp is required for the MCP tools server: pip install sc-tools[mcp]"
    ) from _mcp_err

mcp = FastMCP("sc-tools")

# ---------------------------------------------------------------------------
# Tool: validate_checkpoint
# ---------------------------------------------------------------------------


@mcp.tool()
def validate_checkpoint(uri: str, phase: str, fix: bool = False) -> str:
    """Validate an AnnData checkpoint against the sc_tools metadata contract.

    Parameters
    ----------
    uri
        Local path or URI to the .h5ad checkpoint file.
    phase
        Phase identifier: p1, p2, p3, p35, or p4.
    fix
        If True, attempt auto-fix (e.g. rename obs['batch'] -> obs['raw_data_dir']).

    Returns
    -------
    str
        Validation result summary.
    """
    from sc_tools.validate import validate_file

    cmd = f"validate {phase} {uri}"
    try:
        issues = validate_file(uri, phase=phase, fix=fix)
        result = CLIResult(
            status=Status.success if not issues else Status.error,
            command=cmd,
            data={"issues": issues, "phase": phase, "uri": uri},
            provenance=Provenance(command=cmd),
            message=f"Validation {'passed' if not issues else 'failed'}: {len(issues)} issues found",
        )
        return result.model_dump_json()
    except Exception as exc:
        result = CLIResult(
            status=Status.error,
            command=cmd,
            data={"phase": phase, "uri": uri},
            provenance=Provenance(command=cmd),
            message=f"Validation error: {exc}",
        )
        return result.model_dump_json()


# ---------------------------------------------------------------------------
# Tool: generate_qc_report
# ---------------------------------------------------------------------------


@mcp.tool()
def generate_qc_report(
    report_type: str,
    adata_uri: str,
    figures_dir: str,
    adata_integrated_uri: str = "",
    batch_key: str = "sample",
    viable_min_counts: int = 10,
    viable_min_genes: int = 5,
    soft_min_counts: int = 2,
    neighborhood_k: int = 8,
    tpm_min_counts: int = 100,
    tpm_min_genes: int = 20,
    tpm_min_area_mm2: float = 0.05,
    min_tpm_spots: int = 1000,
) -> str:
    """Generate a date-versioned QC HTML report.

    Parameters
    ----------
    report_type
        One of: pre_filter, post_filter, post_integration, post_celltyping.
    adata_uri
        Path or URI to the input AnnData checkpoint.
    figures_dir
        Output directory for the HTML report.
    adata_integrated_uri
        Optional path to integrated AnnData (for post_integration report).
    batch_key
        Observation column used as batch identifier.
    viable_min_counts
        Minimum total counts per spot to be considered viable (Tier 1, default 10).
        Used for viable region assessment in pre_filter and post_filter reports.
    viable_min_genes
        Minimum genes detected per spot to be considered viable (Tier 1, default 5).
        Used for viable region assessment in pre_filter and post_filter reports.
    soft_min_counts
        Lower bound for contextual (stromal) spot classification (default 2).
        Applied when spatial coordinates are present. Spots with counts between
        soft_min_counts and viable_min_counts whose k-NN neighborhood median
        counts >= viable_min_counts are classified as stromal-contextual (amber).
    neighborhood_k
        Number of spatial nearest neighbors for contextual classification (default 8).
    tpm_min_counts
        Minimum total counts per spot to be considered TPM-worthy (Tier 2, default 100).
    tpm_min_genes
        Minimum genes detected per spot to be considered TPM-worthy (Tier 2, default 20).
    tpm_min_area_mm2
        Minimum TPM-worthy tissue area in mm² for the insufficient_area failure mode
        (default 0.05).
    min_tpm_spots
        Minimum number of TPM-worthy spots to pass QC (default 1000). Samples below
        this threshold hard-fail; samples above this threshold are rescued from
        median-based failures.

    Returns
    -------
    str
        Shell command to run for this report type.
    """
    cmd_parts = [
        "python scripts/run_qc_report.py",
        f"--adata {adata_uri}",
        f"--figures-dir {figures_dir}",
        f"--report {report_type}",
        f"--batch-key {batch_key}",
    ]
    if adata_integrated_uri:
        cmd_parts.append(f"--adata-integrated {adata_integrated_uri}")
    if report_type in ("pre_filter", "post_filter"):
        cmd_parts.append(f"--viable-min-counts {viable_min_counts}")
        cmd_parts.append(f"--viable-min-genes {viable_min_genes}")
        cmd_parts.append(f"--soft-min-counts {soft_min_counts}")
        cmd_parts.append(f"--neighborhood-k {neighborhood_k}")
        cmd_parts.append(f"--tpm-min-counts {tpm_min_counts}")
        cmd_parts.append(f"--tpm-min-genes {tpm_min_genes}")
        cmd_parts.append(f"--tpm-min-area-mm2 {tpm_min_area_mm2}")
        cmd_parts.append(f"--min-tpm-spots {min_tpm_spots}")
    cmd = " \\\n    ".join(cmd_parts)
    return f"Run the following command to generate the QC report:\n\n{cmd}"


# ---------------------------------------------------------------------------
# Tool: score_gene_signatures
# ---------------------------------------------------------------------------


@mcp.tool()
def score_gene_signatures(
    adata_uri: str,
    gene_sets_uri: str,
    output_uri: str,
    method: str = "scanpy",
) -> str:
    """Score gene signatures on an AnnData checkpoint.

    Parameters
    ----------
    adata_uri
        Path or URI to the input AnnData (p3 or p2 checkpoint).
    gene_sets_uri
        Path or URI to the gene signatures JSON file.
    output_uri
        Path or URI for the output scored AnnData (p35 checkpoint).
    method
        Scoring method: scanpy (default), ucell, or ssgsea.

    Returns
    -------
    str
        Shell command to run gene signature scoring.
    """
    cmd = (
        f'python -c "\n'
        f"import anndata as ad\n"
        f"from sc_tools.storage import smart_read_h5ad, smart_write_checkpoint\n"
        f"from sc_tools.tl import score_signature\n"
        f"from sc_tools.tl.gene_sets import load_hallmark, load_msigdb_json, merge_gene_signatures\n"
        f"adata = smart_read_h5ad('{adata_uri}')\n"
        f"import json\n"
        f"with open('{gene_sets_uri}') as f:\n"
        f"    project_sigs = json.load(f)\n"
        f"gene_sets = merge_gene_signatures(project_sigs, load_hallmark())\n"
        f"score_signature(adata, gene_sets, method='{method}')\n"
        f"smart_write_checkpoint(adata, '{output_uri}')\n"
        f"print('Done: {output_uri}')\n"
        f'"'
    )
    return f"Run the following command to score gene signatures:\n\n{cmd}"


# ---------------------------------------------------------------------------
# Tool: run_deconvolution
# ---------------------------------------------------------------------------


@mcp.tool()
def run_deconvolution(
    adata_uri: str,
    reference_uri: str,
    output_uri: str,
    method: str = "cell2location",
    library_key: str = "library_id",
) -> str:
    """Run cell-type deconvolution on a spatial AnnData.

    Parameters
    ----------
    adata_uri
        Path or URI to the spatial AnnData (p35 or p4 checkpoint).
    reference_uri
        Path or URI to the reference AnnData with cell-type labels.
    output_uri
        Path or URI for the output deconvolved AnnData.
    method
        Backend: cell2location (default), tangram, or destvi.
    library_key
        obs column identifying spatial libraries for per-library batching.

    Returns
    -------
    str
        Shell command to run deconvolution.
    """
    cmd = (
        f'python -c "\n'
        f"from sc_tools.storage import smart_read_h5ad, smart_write_checkpoint\n"
        f"from sc_tools.tl import deconvolution\n"
        f"adata = smart_read_h5ad('{adata_uri}')\n"
        f"ref = smart_read_h5ad('{reference_uri}')\n"
        f"deconvolution(adata, ref, method='{method}', library_key='{library_key}')\n"
        f"smart_write_checkpoint(adata, '{output_uri}')\n"
        f"print('Done: {output_uri}')\n"
        f'"'
    )
    return f"Run the following command to run deconvolution:\n\n{cmd}"


# ---------------------------------------------------------------------------
# Tool: generate_sbatch_spaceranger
# ---------------------------------------------------------------------------


@mcp.tool()
def generate_sbatch_spaceranger(
    sample_id: str,
    fastq_dir: str,
    image: str,
    slide: str,
    area: str,
    output_dir: str,
    cluster: str = "brb",
    modality: str = "visium",
) -> str:
    """Generate a Space Ranger SLURM sbatch script.

    Parameters
    ----------
    sample_id
        Sample identifier.
    fastq_dir
        Path to the directory containing FASTQs.
    image
        Path to the brightfield H&E TIFF image.
    slide
        Visium slide serial number.
    area
        Visium capture area identifier (e.g. A1, B1).
    output_dir
        Output directory for Space Ranger results.
    cluster
        HPC cluster: brb or cayuga.
    modality
        visium or visium_hd.

    Returns
    -------
    str
        The generated sbatch script content.
    """
    from sc_tools.ingest.slurm import build_sbatch_header, build_spaceranger_sbatch

    header = build_sbatch_header(
        job_name=f"spaceranger_{sample_id}",
        partition="scu-cpu",
        cpus=16,
        mem_gb=64,
        time_hours=12,
        cluster=cluster,
    )
    body = build_spaceranger_sbatch(
        sample_id=sample_id,
        fastq_dir=fastq_dir,
        image=image,
        slide=slide,
        area=area,
        output_dir=output_dir,
        modality=modality,
    )
    script = header + "\n" + body
    return f"Generated sbatch script for {sample_id}:\n\n{script}"


# ---------------------------------------------------------------------------
# Tool: generate_sbatch_imc
# ---------------------------------------------------------------------------


@mcp.tool()
def generate_sbatch_imc(
    sample_id: str,
    mcd_uri: str,
    panel_csv_uri: str,
    output_dir: str,
    cluster: str = "brb",
) -> str:
    """Generate an IMC pipeline SLURM sbatch script.

    Parameters
    ----------
    sample_id
        Sample identifier.
    mcd_uri
        Path to the .mcd file.
    panel_csv_uri
        Path to the panel CSV file.
    output_dir
        Output directory for the IMC pipeline results.
    cluster
        HPC cluster: brb or cayuga.

    Returns
    -------
    str
        The generated sbatch script content.
    """
    from sc_tools.ingest.slurm import build_imc_sbatch, build_sbatch_header

    header = build_sbatch_header(
        job_name=f"imc_{sample_id}",
        partition="scu-cpu",
        cpus=8,
        mem_gb=32,
        time_hours=4,
        cluster=cluster,
    )
    body = build_imc_sbatch(
        sample_id=sample_id,
        mcd_file=mcd_uri,
        panel_csv=panel_csv_uri,
        output_dir=output_dir,
    )
    script = header + "\n" + body
    return f"Generated IMC sbatch script for {sample_id}:\n\n{script}"


# ---------------------------------------------------------------------------
# Tool: collect_batch_manifests
# ---------------------------------------------------------------------------


@mcp.tool()
def collect_batch_manifests(metadata_dir: str) -> str:
    """Collect per-batch TSV files into all_samples.tsv.

    Parameters
    ----------
    metadata_dir
        Path to the metadata/phase0/ directory containing *_samples.tsv files.

    Returns
    -------
    str
        Summary of collected samples.
    """
    from sc_tools.ingest.config import collect_all_batches

    try:
        df = collect_all_batches(metadata_dir)
        if df.empty:
            return f"No batch TSV files found in {metadata_dir}."
        summary = f"Collected {len(df)} samples from {metadata_dir}/all_samples.tsv.\n"
        if "sample_id" in df.columns:
            summary += "\nSamples:\n"
            for sid in df["sample_id"].tolist():
                summary += f"  - {sid}\n"
        return summary
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: load_sample
# ---------------------------------------------------------------------------


@mcp.tool()
def load_sample(modality: str, sample_id: str, data_uri: str) -> str:
    """Load one sample into AnnData and report its shape.

    Parameters
    ----------
    modality
        One of: visium, visium_hd, visium_hd_cell, xenium, imc.
    sample_id
        Sample identifier.
    data_uri
        Path or URI to the sample data directory or h5ad file.

    Returns
    -------
    str
        Summary of loaded AnnData (shape, obs keys, obsm keys).
    """
    from sc_tools.storage import with_local_copy

    loaders = {
        "visium": "load_visium_sample",
        "visium_hd": "load_visium_hd_sample",
        "visium_hd_cell": "load_visium_hd_cell_sample",
        "xenium": "load_xenium_sample",
        "imc": "load_imc_sample",
    }

    if modality not in loaders:
        return f"ERROR: Unknown modality '{modality}'. Must be one of: {list(loaders)}"

    loader_name = loaders[modality]

    try:
        import sc_tools.ingest as ingest

        loader = getattr(ingest, loader_name)
        with with_local_copy(data_uri) as local:
            adata = loader(local, sample_id)
        return (
            f"Loaded {modality} sample '{sample_id}': "
            f"{adata.n_obs} obs x {adata.n_vars} vars\n"
            f"  obs keys : {list(adata.obs.columns[:10])}\n"
            f"  obsm keys: {list(adata.obsm.keys())}"
        )
    except Exception as exc:
        return f"ERROR loading {sample_id}: {exc}"


# ---------------------------------------------------------------------------
# Tool: inspect_checkpoint
# ---------------------------------------------------------------------------


@mcp.tool()
def inspect_checkpoint(uri: str) -> str:
    """Inspect an AnnData checkpoint structure without loading arrays into memory.

    Uses h5py to read only metadata (shape, column names, keys) so it is fast
    even for large files on HPC scratch.

    Parameters
    ----------
    uri
        Local path or URI to the .h5ad file.

    Returns
    -------
    str
        Summary of obs shape, obs columns, obsm keys, uns keys, and layers.
    """
    from sc_tools.storage import with_local_copy

    try:
        import h5py
    except ImportError:
        return "ERROR: h5py is required for inspect_checkpoint. Install with: pip install h5py"

    try:
        with with_local_copy(uri) as local:
            with h5py.File(local, "r") as f:
                # Shape
                n_obs = f["obs"].attrs.get("_index", None)
                if n_obs is None and "_index" in f["obs"]:
                    n_obs = len(f["obs"]["_index"])
                elif n_obs is None:
                    first_key = next(iter(f["obs"].keys()), None)
                    n_obs = len(f["obs"][first_key]) if first_key else "?"

                n_vars = f["var"].attrs.get("_index", None)
                if n_vars is None and "_index" in f["var"]:
                    n_vars = len(f["var"]["_index"])
                elif n_vars is None:
                    first_key = next(iter(f["var"].keys()), None)
                    n_vars = len(f["var"][first_key]) if first_key else "?"

                # obs columns (skip internal _index)
                obs_cols = [k for k in f["obs"].keys() if not k.startswith("__")]

                # obsm keys
                obsm_keys = list(f["obsm"].keys()) if "obsm" in f else []

                # uns keys (top level only)
                uns_keys = list(f["uns"].keys()) if "uns" in f else []

                # layers
                layer_keys = list(f["layers"].keys()) if "layers" in f else []

                # obsp
                obsp_keys = list(f["obsp"].keys()) if "obsp" in f else []

        lines = [
            f"Checkpoint: {uri}",
            "=" * 40,
            f"  Shape    : {n_obs} obs x {n_vars} vars",
            f"  obs cols : {', '.join(obs_cols[:20])}{'...' if len(obs_cols) > 20 else ''}",
            f"  obsm     : {', '.join(obsm_keys) or 'none'}",
            f"  layers   : {', '.join(layer_keys) or 'none'}",
            f"  obsp     : {', '.join(obsp_keys) or 'none'}",
            f"  uns      : {', '.join(uns_keys[:15])}{'...' if len(uns_keys) > 15 else ''}",
        ]
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_signatures
# ---------------------------------------------------------------------------


@mcp.tool()
def list_signatures(gene_sets_uri: str) -> str:
    """List gene signatures available in a JSON gene sets file.

    Parameters
    ----------
    gene_sets_uri
        Local path or URI to a gene signatures JSON file.
        Accepts flat ``{name: [genes]}`` or nested ``{group: {name: [genes]}}`` format.

    Returns
    -------
    str
        Table of signature names with gene counts, grouped by pathway/category.
    """
    from sc_tools.storage import with_local_copy

    try:
        import json

        with with_local_copy(gene_sets_uri) as local:
            with open(local) as fh:
                data = json.load(fh)

        lines = [f"Gene signatures in: {gene_sets_uri}", "=" * 40]
        total_sigs = 0

        for key, val in data.items():
            if isinstance(val, list):
                # Flat format: {name: [genes]}
                n = len(val)
                lines.append(f"  {key:<40} {n:>4} genes")
                total_sigs += 1
            elif isinstance(val, dict):
                # Nested format: {group: {name: [genes]}}
                lines.append(f"\n  [{key}]")
                for sig_name, genes in val.items():
                    n = len(genes) if isinstance(genes, list) else 0
                    lines.append(f"    {sig_name:<38} {n:>4} genes")
                    total_sigs += 1

        lines.append(f"\n  Total: {total_sigs} signatures")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: generate_sbatch_xenium
# ---------------------------------------------------------------------------


@mcp.tool()
def generate_sbatch_xenium(
    sample_id: str,
    xenium_bundle: str,
    output_dir: str,
    cluster: str = "brb",
    xenium_ranger_path: str = "xeniumranger",
) -> str:
    """Generate a Xenium Ranger SLURM sbatch script.

    Parameters
    ----------
    sample_id
        Sample identifier.
    xenium_bundle
        Path to the raw Xenium bundle directory.
    output_dir
        Output directory for Xenium Ranger results.
    cluster
        HPC cluster: brb or cayuga (used to set log directory label).
    xenium_ranger_path
        Path to the xeniumranger binary.

    Returns
    -------
    str
        The generated sbatch script content.
    """
    from sc_tools.ingest.slurm import build_xenium_sbatch

    script = build_xenium_sbatch(
        sample_id=sample_id,
        xenium_bundle=xenium_bundle,
        output_dir=output_dir,
        xenium_ranger_path=xenium_ranger_path,
        log_dir=f"logs/{cluster}",
    )
    return f"Generated Xenium Ranger sbatch script for {sample_id}:\n\n{script}"


# ---------------------------------------------------------------------------
# Tool: environment_info
# ---------------------------------------------------------------------------


@mcp.tool()
def environment_info() -> str:
    """Report the current runtime environment: package versions and compute backends.

    Returns
    -------
    str
        Formatted summary of Python version, key package versions,
        GPU devices (if any), and availability of rapids/dask.
    """
    from sc_tools.memory.env import get_environment_info

    info = get_environment_info()

    lines = [
        "Runtime Environment",
        "=" * 40,
        f"  Python : {info['python'].split()[0]}",
        "",
        "  Packages:",
    ]
    for pkg, ver in info["packages"].items():
        status = ver if ver is not None else "not installed"
        lines.append(f"    {pkg:<24} {status}")

    lines.append("")
    lines.append("  Compute backends:")

    if info["gpu"]:
        for dev in info["gpu"]:
            lines.append(f"    GPU [{dev['index']}] {dev['name']}  {dev['memory_gb']} GB")
    else:
        lines.append("    GPU              not available")

    lines.append(f"    rapids-singlecell  {'available' if info['rapids'] else 'not available'}")
    lines.append(f"    dask               {'available' if info['dask'] else 'not available'}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Tool: run_full_phase
# ---------------------------------------------------------------------------


@mcp.tool()
def run_full_phase(
    project_path: str,
    phase_slug: str,
    dry_run: bool = False,
) -> dict:
    """Prepare and return the Snakemake command for a full pipeline phase.

    Does not execute the command -- returns it as a string for the agent to
    run via Bash after reviewing.  Validates that the project directory exists
    and that all prerequisite checkpoint files are present before returning.

    Parameters
    ----------
    project_path
        Absolute path to the project directory (must contain a Snakefile).
    phase_slug
        Pipeline phase identifier, e.g. qc_filter, metadata_attach, preprocess,
        scoring, celltype_manual, biology, meta_analysis, ingest_raw, ingest_load.
    dry_run
        If True, return a description of what would be run without changing
        any files.  If False (default), return the ready-to-run command string.

    Returns
    -------
    dict
        Keys: phase, status, command, inputs, outputs, missing, message.
        status is one of: ready, missing_inputs, error.
    """
    from pathlib import Path

    from sc_tools.pipeline import get_phase, get_phase_checkpoint

    # ------------------------------------------------------------------
    # 1. Validate project directory
    # ------------------------------------------------------------------
    proj = Path(project_path)
    if not proj.exists():
        return {
            "phase": phase_slug,
            "status": "error",
            "command": "",
            "inputs": [],
            "outputs": [],
            "missing": [],
            "message": f"Project path does not exist: {project_path}",
        }
    if not proj.is_dir():
        return {
            "phase": phase_slug,
            "status": "error",
            "command": "",
            "inputs": [],
            "outputs": [],
            "missing": [],
            "message": f"Project path is not a directory: {project_path}",
        }

    # ------------------------------------------------------------------
    # 2. Validate phase slug
    # ------------------------------------------------------------------
    try:
        spec = get_phase(phase_slug)
    except KeyError as exc:
        return {
            "phase": phase_slug,
            "status": "error",
            "command": "",
            "inputs": [],
            "outputs": [],
            "missing": [],
            "message": f"Unknown phase slug '{phase_slug}': {exc}",
        }

    # ------------------------------------------------------------------
    # 3. Determine input checkpoints (checkpoints of depends_on phases)
    # ------------------------------------------------------------------
    inputs: list[str] = []
    for dep_slug in spec.depends_on:
        try:
            dep_checkpoint = get_phase_checkpoint(dep_slug)
        except KeyError:
            dep_checkpoint = None
        if dep_checkpoint is not None:
            inputs.append(dep_checkpoint)

    # ------------------------------------------------------------------
    # 4. Check that input checkpoints exist on disk
    # ------------------------------------------------------------------
    missing: list[str] = []
    for rel_path in inputs:
        abs_path = proj / rel_path
        if not abs_path.exists():
            missing.append(rel_path)

    # ------------------------------------------------------------------
    # 5. Determine output checkpoint(s) for this phase
    # ------------------------------------------------------------------
    try:
        phase_checkpoint = get_phase_checkpoint(phase_slug)
        outputs: list[str] = [phase_checkpoint] if phase_checkpoint is not None else []
    except KeyError:
        outputs = []
        # placeholder phase (e.g. ingest_load/{sample_id}) — outputs not resolvable without sample_id

    # ------------------------------------------------------------------
    # 6. Locate Snakefile
    # ------------------------------------------------------------------
    snakefile = proj / "Snakefile"
    snakefile_path = str(snakefile)

    # ------------------------------------------------------------------
    # 7. Build Snakemake command
    # ------------------------------------------------------------------
    command = f"snakemake -d {project_path} -s {snakefile_path} {phase_slug}"
    if dry_run:
        command += " --dry-run"

    # ------------------------------------------------------------------
    # 8. Return structured result
    # ------------------------------------------------------------------
    status = "missing_inputs" if missing else "ready"
    if dry_run:
        message = "Dry run: command not executed."
    elif not missing:
        message = f"Phase {phase_slug} is ready to run."
    else:
        message = f"Missing prerequisite checkpoints: {missing}"

    return {
        "phase": phase_slug,
        "status": status,
        "command": command,
        "inputs": inputs,
        "outputs": outputs,
        "missing": missing,
        "message": message,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mcp.run()
