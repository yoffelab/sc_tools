"""Technology taxonomy + per-phase tracking.

Revision ID: 0002
Revises: 0001
Create Date: 2026-03-06

Changes
-------
* ``projects`` table: add ``domain`` and ``imaging_modality`` columns.
* ``datasets`` table: add ``file_role``, ``validated``, ``n_obs``, ``n_vars`` columns.
* New ``project_phases`` table: per-phase pipeline status with composite PK
  ``(project_id, phase)``.

Domain / imaging_modality vocabulary
--------------------------------------
domain:
    spatial_transcriptomics | spatial_proteomics | imaging | single_cell | bulk

imaging_modality:
    brightfield | fluorescence | multiplexed_fluorescence | probe_based |
    mass_spec_imaging | sequencing_based

file_role:
    primary | supplementary | entry_point | spatialdata | image | metadata
"""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op

revision = "0002"
down_revision = "0001"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ── projects: add domain and imaging_modality ─────────────────────────────
    with op.batch_alter_table("projects") as batch_op:
        batch_op.add_column(sa.Column("domain", sa.String(), nullable=True))
        batch_op.add_column(sa.Column("imaging_modality", sa.String(), nullable=True))

    # ── datasets: add file_role, validated, n_obs, n_vars ─────────────────────
    with op.batch_alter_table("datasets") as batch_op:
        batch_op.add_column(
            sa.Column("file_role", sa.String(), nullable=True, server_default="primary")
        )
        batch_op.add_column(sa.Column("validated", sa.Boolean(), nullable=True, server_default="0"))
        batch_op.add_column(sa.Column("n_obs", sa.Integer(), nullable=True))
        batch_op.add_column(sa.Column("n_vars", sa.Integer(), nullable=True))

    # ── project_phases: new table ─────────────────────────────────────────────
    op.create_table(
        "project_phases",
        sa.Column("project_id", sa.Integer(), nullable=False),
        sa.Column("phase", sa.String(), nullable=False),
        sa.Column("status", sa.String(), nullable=True, server_default="not_started"),
        sa.Column("entry_phase", sa.Boolean(), nullable=True, server_default="0"),
        sa.Column("primary_dataset_id", sa.Integer(), nullable=True),
        sa.Column("n_obs", sa.Integer(), nullable=True),
        sa.Column("n_vars", sa.Integer(), nullable=True),
        sa.Column("n_samples", sa.Integer(), nullable=True),
        sa.Column("notes", sa.Text(), nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
        sa.Column("updated_at", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["project_id"], ["projects.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["primary_dataset_id"], ["datasets.id"]),
        sa.PrimaryKeyConstraint("project_id", "phase"),
    )


def downgrade() -> None:
    op.drop_table("project_phases")

    with op.batch_alter_table("datasets") as batch_op:
        batch_op.drop_column("n_vars")
        batch_op.drop_column("n_obs")
        batch_op.drop_column("validated")
        batch_op.drop_column("file_role")

    with op.batch_alter_table("projects") as batch_op:
        batch_op.drop_column("imaging_modality")
        batch_op.drop_column("domain")
