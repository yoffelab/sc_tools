"""Rename project_data_inventory → data_project_map, add patient_data_map, drop data_phase.patient_id FK.

Revision ID: 0010
Revises: 0009
Create Date: 2026-03-13

Changes
-------
* Rename ``project_data_inventory`` → ``data_project_map``
* Create ``patient_data_map`` join table (many-to-many patients ↔ data)
* Migrate existing ``data_phase.patient_id`` links into ``patient_data_map``
* Drop ``patient_id`` column from ``data_phase``
"""

from __future__ import annotations

revision = "0010"
down_revision = "0009"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # 1. Rename join table
    op.rename_table("project_data_inventory", "data_project_map")

    # 2. Create patient_data_map
    op.create_table(
        "patient_data_map",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column(
            "patient_id",
            sa.Integer,
            sa.ForeignKey("patients.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "data_id",
            sa.Integer,
            sa.ForeignKey("data_phase.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("role", sa.String, server_default="source"),
        sa.Column("notes", sa.Text),
    )

    # 3. Migrate existing patient_id links from data_phase into patient_data_map
    op.execute(
        """
        INSERT INTO patient_data_map (patient_id, data_id, role)
        SELECT patient_id, id, 'source'
        FROM data_phase
        WHERE patient_id IS NOT NULL
        """
    )

    # 4. Drop patient_id column from data_phase
    with op.batch_alter_table("data_phase") as batch_op:
        batch_op.drop_column("patient_id")


def downgrade() -> None:
    # Re-add patient_id column
    with op.batch_alter_table("data_phase") as batch_op:
        batch_op.add_column(
            sa.Column("patient_id", sa.Integer, sa.ForeignKey("patients.id"), nullable=True)
        )

    # Migrate links back
    op.execute(
        """
        UPDATE data_phase SET patient_id = pdm.patient_id
        FROM patient_data_map pdm
        WHERE data_phase.id = pdm.data_id
        """
    )

    op.drop_table("patient_data_map")
    op.rename_table("data_project_map", "project_data_inventory")
