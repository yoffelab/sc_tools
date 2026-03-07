"""Alembic environment configuration for sc_tools registry.

This file is loaded by ``alembic upgrade`` / ``alembic revision``.
It connects to the database and runs migrations.

Configure ``alembic.ini`` at repo root::

    [alembic]
    script_location = sc_tools/migrations
    sqlalchemy.url = sqlite:///%(here)s/.sc_tools/registry.db

Override the URL at runtime::

    SC_TOOLS_REGISTRY_URL=postgresql://... alembic upgrade head
"""

from __future__ import annotations

import os
from logging.config import fileConfig

from alembic import context

# ---------------------------------------------------------------------------
# Import the registry Base so Alembic can detect model changes
# ---------------------------------------------------------------------------
# We build models on the fly in registry.py, so we reproduce the Base here.
try:
    from sqlalchemy.orm import DeclarativeBase

    class Base(DeclarativeBase):
        pass

    # Re-import model builder so Alembic sees the table definitions
    from sc_tools.registry import _build_models

    _build_models(Base)
    target_metadata = Base.metadata
except Exception:
    target_metadata = None  # type: ignore[assignment]

# ---------------------------------------------------------------------------
# Alembic Config object
# ---------------------------------------------------------------------------
config = context.config

if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Allow SC_TOOLS_REGISTRY_URL env var to override alembic.ini
db_url = os.environ.get("SC_TOOLS_REGISTRY_URL")
if db_url:
    config.set_main_option("sqlalchemy.url", db_url)


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode (no live DB connection)."""
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )
    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations with a live DB connection."""
    from sqlalchemy import engine_from_config, pool

    connectable = engine_from_config(
        config.get_section(config.config_ini_section, {}),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )
    with connectable.connect() as connection:
        context.configure(connection=connection, target_metadata=target_metadata)
        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
