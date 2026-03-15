#!/usr/bin/env python3
"""One-shot migration: copy all data from SQLite registry to PostgreSQL (Supabase).

Usage
-----
    # Dry-run (print row counts, don't write):
    python scripts/migrate_sqlite_to_pg.py --dry-run

    # Migrate to local Supabase:
    python scripts/migrate_sqlite_to_pg.py \
        --target "postgresql://postgres:postgres@127.0.0.1:54322/postgres"

    # Migrate to cloud Supabase:
    python scripts/migrate_sqlite_to_pg.py \
        --target "postgresql://postgres.<ref>:<pw>@aws-0-us-east-1.pooler.supabase.com:6543/postgres"

Assumes Alembic migrations have already been run on the target database.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from sqlalchemy import MetaData, create_engine, text
from sqlalchemy.orm import Session

# SQLite source — default location
DEFAULT_SQLITE = f"sqlite:///{Path.home() / '.sc_tools' / 'registry.db'}"

# Tables in dependency order (parents before children)
TABLE_ORDER = [
    "projects",
    "data_sources",
    "project_data_sources",
    "subjects",
    "samples",
    "subject_project_links",
    "project_phases",
    "datasets",
    "bio_data",
    "bio_images",
    "rnaseq_data",
    "spatial_seq_data",
    "epigenomics_data",
    "genome_seq_data",
    "slurm_jobs",
    "agent_tasks",
]

# Alembic's own table — skip during migration
SKIP_TABLES = {"alembic_version"}


def migrate(source_url: str, target_url: str, *, dry_run: bool = False) -> None:
    src_engine = create_engine(source_url)
    tgt_engine = create_engine(target_url)

    # Reflect source schema
    src_meta = MetaData()
    src_meta.reflect(bind=src_engine)

    src_tables = {t.name for t in src_meta.sorted_tables} - SKIP_TABLES

    # Determine migration order: use TABLE_ORDER for known tables, append any extras
    ordered = [t for t in TABLE_ORDER if t in src_tables]
    extras = sorted(src_tables - set(ordered))
    ordered.extend(extras)

    print(f"Source: {source_url}")
    print(f"Target: {target_url}")
    print(f"Tables to migrate: {len(ordered)}")
    print()

    total_rows = 0

    with Session(src_engine) as src_session, Session(tgt_engine) as tgt_session:
        if not dry_run:
            # Disable FK checks during bulk load (superuser required — Supabase local has it)
            tgt_session.execute(text("SET session_replication_role = 'replica'"))

        for table_name in ordered:
            table = src_meta.tables[table_name]
            rows = src_session.execute(table.select()).mappings().all()
            count = len(rows)
            total_rows += count
            print(f"  {table_name}: {count} rows", end="")

            if dry_run or count == 0:
                print(" (skip)" if dry_run else "")
                continue

            # Reflect the target table to get correct column types
            tgt_meta = MetaData()
            tgt_meta.reflect(bind=tgt_engine, only=[table_name])
            tgt_table = tgt_meta.tables[table_name]

            # Clear existing rows in target (idempotent re-runs)
            tgt_session.execute(tgt_table.delete())

            # Insert in batches of 500
            batch_size = 500
            for i in range(0, count, batch_size):
                batch = [dict(r) for r in rows[i : i + batch_size]]
                tgt_session.execute(tgt_table.insert(), batch)

            print(" ✓")

        if not dry_run:
            # Re-enable FK checks before commit
            tgt_session.execute(text("SET session_replication_role = 'origin'"))
            tgt_session.commit()
            print("\nData committed.")

            # Reset PostgreSQL sequences to max(id)+1
            print("\nResetting sequences...")
            for table_name in ordered:
                table = src_meta.tables[table_name]
                # Only reset for tables with an integer 'id' primary key
                pk_cols = [c for c in table.primary_key.columns]
                if len(pk_cols) == 1 and pk_cols[0].name == "id":
                    seq_name = f"{table_name}_id_seq"
                    try:
                        result = tgt_session.execute(
                            text(f"SELECT MAX(id) FROM {table_name}")
                        ).scalar()
                        if result is not None:
                            next_val = result + 1
                            tgt_session.execute(
                                text(f"SELECT setval('{seq_name}', {next_val}, false)")
                            )
                            print(f"  {seq_name} → {next_val}")
                    except Exception as e:
                        print(f"  {seq_name}: skipped ({e})")
                        tgt_session.rollback()

            tgt_session.commit()

    print(f"\nDone. {total_rows} total rows {'would be ' if dry_run else ''}migrated.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Migrate sc_tools registry from SQLite to PostgreSQL")
    parser.add_argument(
        "--source",
        default=DEFAULT_SQLITE,
        help=f"SQLite source URL (default: {DEFAULT_SQLITE})",
    )
    parser.add_argument(
        "--target",
        default="postgresql://postgres:postgres@127.0.0.1:54322/postgres",
        help="PostgreSQL target URL (default: local Supabase)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print row counts without writing",
    )
    args = parser.parse_args()

    if "sqlite" not in args.source:
        print("ERROR: --source should be a SQLite URL", file=sys.stderr)
        sys.exit(1)
    if "postgresql" not in args.target and "postgres" not in args.target:
        print("ERROR: --target should be a PostgreSQL URL", file=sys.stderr)
        sys.exit(1)

    migrate(args.source, args.target, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
