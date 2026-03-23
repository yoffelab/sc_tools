"""CLIResult Pydantic envelope for structured CLI and MCP output.

Every ``sct`` command returns a :class:`CLIResult` instance.  The same model
is used by the MCP server, providing dual serialization:

* **CLI:** ``result.model_dump_json()`` -- JSON string to stdout.
* **MCP:** ``result.model_dump(mode="json")`` -- JSON-safe dict for tool returns.

Fields follow decisions D-06 through D-11 from the Phase 2 context.
"""

from __future__ import annotations

import importlib.metadata
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field


def _get_version() -> str:
    """Return sc_tools package version, or ``'unknown'`` if not installed."""
    try:
        return importlib.metadata.version("sci-sc-tools")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


class Status(str, Enum):
    """Command execution status (D-07)."""

    success = "success"
    error = "error"
    partial = "partial"
    skipped = "skipped"


class ErrorInfo(BaseModel):
    """Structured error information with taxonomy (D-13, D-14).

    ``category`` is one of ``"retryable"``, ``"fixable"``, or ``"fatal"``.
    ``suggestion`` provides an actionable fix the agent can act on directly.
    """

    category: str  # "retryable" | "fixable" | "fatal"
    suggestion: str
    details: str | None = None


class InputFile(BaseModel):
    """Input file record for provenance tracking (D-08).

    Records the path, type, SHA256 checksum, and size of each input file
    consumed by a CLI command.
    """

    path: str
    path_type: str = "relative"
    sha256: str
    size_bytes: int


class ProvenanceRecord(BaseModel):
    """Full provenance sidecar content (D-06).

    Used for .provenance.json sidecar files and adata.uns embedding.
    Separate from the minimal Provenance class used in CLIResult for
    backwards compatibility.
    """

    command: str
    params: dict[str, Any] = Field(default_factory=dict)
    inputs: list[InputFile] = Field(default_factory=list)
    sc_tools_version: str = Field(default_factory=_get_version)
    timestamp: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),  # noqa: UP017
    )
    runtime_s: float | None = None
    peak_memory_mb: float | None = None


class Provenance(BaseModel):
    """Minimal provenance metadata (D-10).

    Kept for backwards compatibility in CLIResult. Full provenance uses
    ProvenanceRecord for sidecars.
    """

    command: str
    timestamp: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),  # noqa: UP017
    )
    sc_tools_version: str = Field(default_factory=_get_version)


class CLIResult(BaseModel):
    """Structured envelope returned by every ``sct`` command (D-06).

    Attributes
    ----------
    status : Status
        Execution outcome.
    command : str
        The CLI command that produced this result.
    data : dict
        Computed results (metrics, status tables).
    artifacts : list[str]
        File paths created (h5ad, HTML reports, plots).
    provenance : Provenance
        Command metadata (timestamp, version).
    message : str
        Human-readable summary of what happened.
    error : ErrorInfo | None
        Structured error info when ``status`` is ``error``.
    partial_failures : list[dict] | None
        Per-item failure details when ``status`` is ``partial`` (D-08).
        Structured so agents can identify and retry only failed items.
    """

    model_config = ConfigDict(use_enum_values=True)

    status: Status
    command: str
    data: dict[str, Any] = Field(default_factory=dict)
    artifacts: list[str] = Field(default_factory=list)
    provenance: Provenance
    message: str = ""
    error: ErrorInfo | None = None
    partial_failures: list[dict[str, Any]] | None = None
