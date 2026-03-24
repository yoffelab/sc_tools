"""Estimate command for pre-execution memory/runtime projection (MEM-02).

Provides ``sct estimate <command_name> <file>`` to estimate peak memory
and runtime before running expensive operations.
"""

from __future__ import annotations

import typer


def register_estimate(app: typer.Typer) -> None:
    """Register the estimate command on the given Typer app."""
    from sc_tools.cli import cli_handler

    @app.command("estimate")
    @cli_handler
    def estimate_command(
        command_name: str = typer.Argument(
            ..., help="Command to estimate (e.g., preprocess, qc, de)"
        ),
        file: str = typer.Argument(..., help="Input h5ad file path"),
    ) -> None:
        """Estimate peak memory and runtime for a command before execution."""
        from pathlib import Path

        from sc_tools.errors import SCToolsUserError
        from sc_tools.io.estimate import METHOD_MULTIPLIERS, estimate_from_h5
        from sc_tools.models.result import CLIResult, Provenance, Status

        file_path = Path(file)
        if not file_path.exists():
            raise SCToolsUserError(
                f"File not found: {file}",
                suggestion="Check the file path",
            )

        est = estimate_from_h5(str(file_path))

        # Apply method-specific multiplier
        method_key = command_name.split()[0]  # Handle "preprocess run" -> "preprocess"
        multiplier = METHOD_MULTIPLIERS.get(method_key, METHOD_MULTIPLIERS["default"])
        est["method_multiplier"] = multiplier
        est["estimated_peak_mb"] = est["base_mb"] * multiplier

        # Runtime estimate (rough): cells per second varies by operation
        throughput = {
            "preprocess": 5000,
            "de": 8000,
            "qc": 20000,
            "benchmark": 10000,
        }
        cells_per_sec = throughput.get(method_key, 10000)
        est["estimated_runtime_s"] = round(est["n_obs"] / cells_per_sec, 1)
        est["runtime_confidence"] = "low"
        est["memory_confidence"] = "high"

        return CLIResult(
            status=Status.success,
            command=f"estimate {command_name}",
            data=est,
            provenance=Provenance(command=f"estimate {command_name}"),
            message=(
                f"Estimated peak memory: {est['estimated_peak_mb']:.0f}MB "
                f"for {est['n_obs']} cells x {est['n_vars']} genes"
            ),
        )
