"""Tests for Alembic migration chain integrity.

Verifies that ``alembic upgrade head`` produces the expected schema and that
the ORM (``create_all``) and migration chain stay in sync.
"""

from __future__ import annotations

import pytest

sa = pytest.importorskip("sqlalchemy")
pytest.importorskip("alembic")

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

EXPECTED_TABLES = {
    "projects",
    "datasets",
    "slurm_jobs",
    "agent_tasks",
    "project_phases",
    "data_sources",
    "project_data_sources",
    "subjects",
    "samples",
    "subject_project_links",
    "bio_data",
    "bio_images",
    "rnaseq_data",
    "spatial_seq_data",
    "epigenomics_data",
    "genome_seq_data",
}


def _run_alembic_upgrade(db_url: str) -> None:
    """Run ``alembic upgrade head`` against *db_url* programmatically."""
    import os

    from alembic import command
    from alembic.config import Config

    migrations_dir = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "migrations",
    )

    cfg = Config()
    cfg.set_main_option("script_location", migrations_dir)
    cfg.set_main_option("sqlalchemy.url", db_url)
    command.upgrade(cfg, "head")


def _get_table_names(db_url: str) -> set[str]:
    """Return the set of table names in the database at *db_url*."""
    engine = sa.create_engine(db_url)
    inspector = sa.inspect(engine)
    tables = set(inspector.get_table_names())
    engine.dispose()
    return tables


def _get_column_names(db_url: str, table: str) -> set[str]:
    """Return column names for *table* in the database at *db_url*."""
    engine = sa.create_engine(db_url)
    inspector = sa.inspect(engine)
    columns = {col["name"] for col in inspector.get_columns(table)}
    engine.dispose()
    return columns


# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------


def test_upgrade_head_clean_db(tmp_path):
    """Alembic upgrade head on a fresh DB creates all expected tables."""
    db_path = tmp_path / "test_upgrade.db"
    db_url = f"sqlite:///{db_path}"

    _run_alembic_upgrade(db_url)

    tables = _get_table_names(db_url)
    # alembic_version is always created by Alembic itself
    missing = EXPECTED_TABLES - tables
    assert not missing, f"Missing tables after upgrade head: {missing}"


def test_modality_column_exists(tmp_path):
    """After upgrade head, bio_data table contains a ``modality`` column."""
    db_path = tmp_path / "test_modality.db"
    db_url = f"sqlite:///{db_path}"

    _run_alembic_upgrade(db_url)

    columns = _get_column_names(db_url, "bio_data")
    assert "modality" in columns, f"bio_data missing 'modality' column; has: {sorted(columns)}"


def test_orm_matches_migration(tmp_path):
    """ORM create_all and alembic upgrade head produce the same tables."""
    from sqlalchemy.orm import DeclarativeBase

    from sc_tools.registry import _build_models

    # --- DB via Alembic migrations ---
    alembic_db = tmp_path / "alembic.db"
    alembic_url = f"sqlite:///{alembic_db}"
    _run_alembic_upgrade(alembic_url)
    alembic_tables = _get_table_names(alembic_url) - {"alembic_version"}

    # --- DB via ORM create_all ---
    orm_db = tmp_path / "orm.db"
    orm_url = f"sqlite:///{orm_db}"

    class Base(DeclarativeBase):
        pass

    _build_models(Base)
    engine = sa.create_engine(orm_url)
    Base.metadata.create_all(engine)
    orm_tables = _get_table_names(orm_url)
    engine.dispose()

    assert alembic_tables == orm_tables, (
        f"Table mismatch.\n"
        f"  Only in migrations: {alembic_tables - orm_tables}\n"
        f"  Only in ORM:        {orm_tables - alembic_tables}"
    )
