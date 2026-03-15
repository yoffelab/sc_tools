from __future__ import annotations

revision = "0007"
down_revision = "0006"
branch_labels = None
depends_on = None

"""Consolidate datasets -> bio_data and drop agent_tasks.

Revision ID: 0007
Revises: 0006
Create Date: 2026-03-13

Changes
-------
* Change ``project_phases.primary_dataset_id`` FK from ``datasets.id``
  to ``bio_data.id``. The actual ID values are already identical since
  datasets and bio_data have matching rows (1:1 migration).
* Drop the ``agent_tasks`` table (0 rows, unused).
* The ``datasets`` table is left in place (deprecated) but code no longer
  writes to it.
"""

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # 1. Re-point project_phases.primary_dataset_id FK: datasets -> bio_data
    op.drop_constraint(
        "project_phases_primary_dataset_id_fkey",
        "project_phases",
        type_="foreignkey",
    )
    op.create_foreign_key(
        "project_phases_primary_dataset_id_fkey",
        "project_phases",
        "bio_data",
        ["primary_dataset_id"],
        ["id"],
    )

    # 2. Drop agent_tasks table
    op.drop_table("agent_tasks")


def downgrade() -> None:
    # 1. Recreate agent_tasks table
    op.create_table(
        "agent_tasks",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("task_type", sa.String(), nullable=True),
        sa.Column("project_id", sa.Integer(), sa.ForeignKey("projects.id"), nullable=False),
        sa.Column("status", sa.String(), server_default="queued", nullable=True),
        sa.Column("inputs_json", sa.Text(), nullable=True),
        sa.Column("outputs_json", sa.Text(), nullable=True),
        sa.Column("error", sa.Text(), nullable=True),
        sa.Column("started_at", sa.String(), nullable=True),
        sa.Column("finished_at", sa.String(), nullable=True),
    )

    # 2. Revert project_phases FK back to datasets
    op.drop_constraint(
        "project_phases_primary_dataset_id_fkey",
        "project_phases",
        type_="foreignkey",
    )
    op.create_foreign_key(
        "project_phases_primary_dataset_id_fkey",
        "project_phases",
        "datasets",
        ["primary_dataset_id"],
        ["id"],
    )
