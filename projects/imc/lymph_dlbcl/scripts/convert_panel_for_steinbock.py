"""Convert sc_tools panel CSV to steinbock-compatible format.

sc_tools panel CSV columns:
    channel      -- raw channel name matching txt file header, e.g. CD3(Er170Di)
    Target       -- human-readable protein name, e.g. CD3
    Metal_Tag    -- isotope code, e.g. Er170
    full         -- 1/0 include in intensity measurement (steinbock keep)
    ilastik      -- 1/0 include in ilastik/segmentation features
    raw_channel_name -- same as channel (explicit copy; used for validation)
    notes        -- free text

steinbock panel CSV columns required:
    channel  -- raw channel name (exact match to txt file header)
    name     -- human-readable name shown in outputs
    keep     -- 1/0 include in intensity measurement
    ilastik  -- 1/0 include in ilastik features

Rows with empty channel/raw_channel_name (e.g. markers not found in panel)
are dropped with a warning.

Usage
-----
    python scripts/convert_panel_for_steinbock.py \\
        metadata/panel_immune_t2.csv \\
        data/processed/DLC0002/panel.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def convert_panel(src: Path, dst: Path) -> None:
    """Read sc_tools panel CSV and write steinbock panel CSV."""
    df = pd.read_csv(src)

    required = {"channel", "full", "ilastik"}
    missing = required - set(df.columns)
    if missing:
        print(
            f"ERROR: panel CSV {src} missing columns: {missing}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Use raw_channel_name as the authoritative channel key when present
    if "raw_channel_name" in df.columns:
        channel_col = df["raw_channel_name"].fillna(df["channel"])
    else:
        channel_col = df["channel"]

    # Human-readable name: prefer Target, fall back to channel
    if "Target" in df.columns:
        name_col = df["Target"].fillna(channel_col)
    else:
        name_col = channel_col

    out = pd.DataFrame(
        {
            "channel": channel_col,
            "name": name_col,
            "keep": df["full"].fillna(0).astype(int),
            "ilastik": df["ilastik"].fillna(0).astype(int),
        }
    )

    # Drop rows with empty channel (markers not found in this panel version)
    n_before = len(out)
    out = out[out["channel"].astype(str).str.strip() != ""]
    n_dropped = n_before - len(out)
    if n_dropped > 0:
        print(
            f"WARNING: dropped {n_dropped} rows with empty channel name from {src.name}",
            file=sys.stderr,
        )

    dst.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dst, index=False)
    print(f"Wrote steinbock panel CSV: {dst} ({len(out)} channels)")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("src", help="sc_tools panel CSV (input)")
    parser.add_argument("dst", help="steinbock panel CSV (output)")
    args = parser.parse_args(argv)

    convert_panel(Path(args.src), Path(args.dst))


if __name__ == "__main__":
    main()
