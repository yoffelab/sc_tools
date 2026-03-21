"""sc_tools exception hierarchy with error taxonomy metadata.

Each exception class carries ``category``, ``suggestion``, and ``exit_code``
class attributes that the CLI boundary handler uses to build structured error
responses.

Exit code mapping (per D-15, D-16):
    0 = success (includes partial)
    1 = user error (bad args, missing file)
    2 = data error (file found but wrong format, validation failed)
    3 = runtime error (OOM, computation failure) or fatal
"""

from __future__ import annotations


class SCToolsError(Exception):
    """Base exception for sc_tools."""

    category: str = "fatal"
    suggestion: str = ""
    exit_code: int = 3

    def __init__(
        self,
        message: str,
        *,
        suggestion: str = "",
        details: str | None = None,
    ):
        super().__init__(message)
        self.suggestion = suggestion or self.__class__.suggestion
        self.details = details


class SCToolsUserError(SCToolsError):
    """User-correctable error (bad args, missing file). Exit code 1."""

    category = "fixable"
    exit_code = 1


class SCToolsDataError(SCToolsError):
    """Data validation error (file exists but wrong format). Exit code 2."""

    category = "fixable"
    exit_code = 2


class SCToolsRuntimeError(SCToolsError):
    """Runtime error (OOM, computation failure). Exit code 3."""

    category = "retryable"
    exit_code = 3


class SCToolsFatalError(SCToolsError):
    """Unrecoverable error. Exit code 3."""

    category = "fatal"
    exit_code = 3
