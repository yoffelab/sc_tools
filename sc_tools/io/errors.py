"""Error types for the IO gateway module."""

from __future__ import annotations


class SCToolsRuntimeError(Exception):
    """Raised when a runtime safety check fails (e.g., memory guard).

    Attributes
    ----------
    category : str
        Error taxonomy: "fatal" means the operation cannot proceed.
    exit_code : int
        Suggested CLI exit code (3 = runtime error).
    """

    category: str = "fatal"
    exit_code: int = 3
