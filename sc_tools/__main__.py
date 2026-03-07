"""CLI entry point for sc_tools.

Usage
-----
::

    python -m sc_tools registry status
"""

from __future__ import annotations

import sys


def main() -> None:
    args = sys.argv[1:]

    if not args:
        print("Usage: python -m sc_tools <command> [args]")
        print("")
        print("Commands:")
        print("  registry status  -- show registry status summary")
        sys.exit(1)

    if args[0] == "registry":
        subcmd = args[1] if len(args) > 1 else "status"
        if subcmd == "status":
            from sc_tools.registry import _cli_status

            _cli_status()
        else:
            print(f"Unknown registry sub-command: {subcmd}")
            sys.exit(1)
    else:
        print(f"Unknown command: {args[0]}")
        sys.exit(1)


if __name__ == "__main__":
    main()
