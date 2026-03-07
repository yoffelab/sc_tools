#!/usr/bin/env python
"""CLI wrapper for sc_tools checkpoint validation.

Usage (preferred — semantic slug names, new nomenclature):
    python scripts/validate_checkpoint.py results/adata.raw.h5ad --phase qc_filter
    python scripts/validate_checkpoint.py results/adata.annotated.h5ad --phase metadata_attach
    python scripts/validate_checkpoint.py results/adata.normalized.h5ad --phase preprocess
    python scripts/validate_checkpoint.py results/adata.scored.h5ad --phase scoring
    python scripts/validate_checkpoint.py results/adata.celltyped.h5ad --phase celltype_manual

Usage (with flags):
    python scripts/validate_checkpoint.py results/adata.raw.h5ad --phase qc_filter --fix
    python scripts/validate_checkpoint.py results/adata.raw.h5ad --phase qc_filter --warn-only

Deprecated (old nomenclature — still accepted, emits DeprecationWarning):
    python scripts/validate_checkpoint.py results/adata.raw.p1.h5ad --phase p1
    python scripts/validate_checkpoint.py results/adata.annotated.p2.h5ad --phase p2

Designed for Snakemake shell rules. Exit code 1 on failure unless --warn-only.
"""

from __future__ import annotations

import argparse
import sys

# Semantic slug choices (new nomenclature — preferred)
_SLUG_CHOICES = ["qc_filter", "metadata_attach", "preprocess", "scoring", "celltype_manual"]

# Legacy p-code choices (old nomenclature — deprecated, still accepted)
_LEGACY_CHOICES = ["p1", "p2", "p3", "p35", "p4"]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate an AnnData checkpoint against phase requirements. "
            "Use semantic slug names (new nomenclature). "
            "Legacy p-codes are deprecated (old nomenclature) but still accepted."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Phase identifier reference:\n"
            "  Preferred (new nomenclature):          Legacy (old nomenclature, DEPRECATED):\n"
            "    qc_filter       Phase 1                p1\n"
            "    metadata_attach Phase 2                p2\n"
            "    preprocess      Phase 3                p3\n"
            "    scoring         Phase 3.5b             p35\n"
            "    celltype_manual Phase 4                p4\n"
        ),
    )
    parser.add_argument("path", help="Path to .h5ad file")
    parser.add_argument(
        "--phase",
        required=True,
        choices=_SLUG_CHOICES + _LEGACY_CHOICES,
        metavar="PHASE",
        help=(
            "Phase to validate against. "
            "Preferred slugs (new nomenclature): " + ", ".join(_SLUG_CHOICES) + ". "
            "Deprecated p-codes (old nomenclature): " + ", ".join(_LEGACY_CHOICES) + "."
        ),
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

    # Warn if legacy p-code was passed at the CLI level (validate.py will also
    # emit a DeprecationWarning via warnings.warn, but a clear stderr message
    # is more visible in Snakemake logs).
    if args.phase in _LEGACY_CHOICES:
        print(
            f"NOTE: --phase {args.phase} uses old nomenclature (deprecated). "
            f"Update to the semantic slug equivalent.",
            file=sys.stderr,
        )

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
