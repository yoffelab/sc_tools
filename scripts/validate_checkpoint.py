#!/usr/bin/env python
"""CLI wrapper for sc_tools checkpoint validation.

Usage:
    python scripts/validate_checkpoint.py results/adata.raw.p1.h5ad --phase p1
    python scripts/validate_checkpoint.py results/adata.raw.p1.h5ad --phase p1 --fix
    python scripts/validate_checkpoint.py results/adata.raw.p1.h5ad --phase p1 --warn-only

Designed for Snakemake shell rules. Exit code 1 on failure unless --warn-only.
"""

from __future__ import annotations

import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate an AnnData checkpoint against phase requirements."
    )
    parser.add_argument("path", help="Path to .h5ad file")
    parser.add_argument(
        "--phase",
        required=True,
        choices=["p1", "p2", "p3", "p35", "p4"],
        help="Pipeline phase to validate against",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Attempt auto-fixes (e.g., rename batch -> raw_data_dir)",
    )
    parser.add_argument(
        "--warn-only",
        action="store_true",
        help="Print warnings but exit 0 even on failure",
    )
    args = parser.parse_args()

    from sc_tools.validate import CheckpointValidationError, validate_file

    try:
        issues = validate_file(
            args.path,
            phase=args.phase,
            strict=not args.warn_only,
            fix=args.fix,
        )
    except CheckpointValidationError as e:
        print(f"FAIL: {e}", file=sys.stderr)
        return 1

    if issues:
        for issue in issues:
            print(f"WARNING: {issue}", file=sys.stderr)
        print(
            f"Checkpoint {args.phase} has {len(issues)} warning(s) (--warn-only mode)",
            file=sys.stderr,
        )
    else:
        print(f"OK: {args.path} passes {args.phase} validation")

    return 0


if __name__ == "__main__":
    sys.exit(main())
