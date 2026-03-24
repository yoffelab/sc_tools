"""Drop legacy tables after four-layer schema migration.

Revision ID: 0014
Revises: 0013
Create Date: 2026-03-16

Drops the old tables that were superseded by migration 0013:
* ``data_processing_phase`` -- replaced by inventory_items + project_phases
* ``data_project_map`` -- replaced by project_datasets
* ``data_inventory`` -- replaced by data_sources

Also removes the transitional ``patient_data_map.data_id`` column.
"""

from __future__ import annotations

revision = "0014"
down_revision = "0013"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # 1. Drop the old data_id column from patient_data_map (kept during burn-in)
    with op.batch_alter_table("patient_data_map") as batch_op:
        batch_op.drop_column("data_id")

    # 2. Drop legacy tables in dependency order
    op.drop_table("data_project_map")
    op.drop_table("data_processing_phase")
    op.drop_table("data_inventory")


def downgrade() -> None:
    # Recreate the old tables (empty -- data migration is in 0013 downgrade)

    op.create_table(
        "data_inventory",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("name", sa.String(), unique=True, nullable=False),
        sa.Column("description", sa.Text()),
        sa.Column("uri", sa.Text(), nullable=False),
        sa.Column("platform", sa.String()),
        sa.Column("domain", sa.String()),
        sa.Column("imaging_modality", sa.String()),
        sa.Column("source_type", sa.String()),
        sa.Column("organism", sa.String()),
        sa.Column("tissue", sa.String()),
        sa.Column("disease", sa.String()),
        sa.Column("n_samples", sa.Integer()),
        sa.Column("n_cells", sa.Integer()),
        sa.Column("publication", sa.String()),
        sa.Column("access_notes", sa.Text()),
        sa.Column("status", sa.String(), server_default="available"),
        sa.Column("created_at", sa.String()),
    )

    op.create_table(
        "data_processing_phase",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("phase", sa.String()),
        sa.Column("status", sa.String(), server_default="pending"),
        sa.Column("qc_status", sa.String(), server_default="pending"),
        sa.Column("uri", sa.Text(), nullable=False),
        sa.Column("format", sa.String()),
        sa.Column("platform", sa.String()),
        sa.Column("category", sa.String()),
        sa.Column("file_role", sa.String(), server_default="primary"),
        sa.Column("n_obs", sa.Integer()),
        sa.Column("n_vars", sa.Integer()),
        sa.Column("size_mb", sa.Float()),
        sa.Column("metadata", sa.JSON(), server_default="{}"),
        sa.Column("created_at", sa.String()),
        sa.Column("updated_at", sa.String()),
    )

    op.create_table(
        "data_project_map",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "data_source_id",
            sa.Integer(),
            sa.ForeignKey("data_inventory.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("role", sa.String(), server_default="input"),
        sa.Column("notes", sa.Text()),
    )

    # Re-add data_id column to patient_data_map
    with op.batch_alter_table("patient_data_map") as batch_op:
        batch_op.add_column(sa.Column("data_id", sa.Integer(), nullable=True))
