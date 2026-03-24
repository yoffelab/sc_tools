"""Add samples table, link inventory_items to samples, drop patient_data_map.

Revision ID: 0017
Revises: 0016
Create Date: 2026-03-17

New tables
----------
* ``samples`` -- biological specimens collected from patients

Schema changes
--------------
* ``inventory_items`` gains ``sample_id`` FK to ``samples.id`` (nullable)
* ``patient_data_map`` is dropped (replaced by inventory_items -> samples -> patients)
"""

from __future__ import annotations

revision = "0017"
down_revision = "0016"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # ---- 1. Create samples table ----
    op.create_table(
        "samples",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "patient_id",
            sa.Integer(),
            sa.ForeignKey("patients.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column("sample_id", sa.String(), unique=True, nullable=False),
        sa.Column("tissue", sa.String(), nullable=True),
        sa.Column("collection_date", sa.String(), nullable=True),
        sa.Column("metadata", sa.JSON(), nullable=True, default={}),
        sa.Column("created_at", sa.String(), nullable=True),
    )
    op.create_index("ix_samples_patient_id", "samples", ["patient_id"])

    # ---- 2. Add sample_id FK to inventory_items ----
    op.add_column(
        "inventory_items",
        sa.Column(
            "sample_id",
            sa.Integer(),
            sa.ForeignKey("samples.id", ondelete="SET NULL"),
            nullable=True,
        ),
    )
    op.create_index("ix_inventory_items_sample_id", "inventory_items", ["sample_id"])

    # ---- 3. Migrate patient_data_map data into samples ----
    # For each (patient_id, inventory_id) in patient_data_map, create a sample
    # and link the inventory item to it.
    conn = op.get_bind()
    rows = conn.execute(
        sa.text(
            "SELECT pdm.id, pdm.patient_id, pdm.inventory_id, p.patient_id AS patient_str "
            "FROM patient_data_map pdm "
            "JOIN patients p ON p.id = pdm.patient_id "
            "WHERE pdm.inventory_id IS NOT NULL"
        )
    ).fetchall()

    from datetime import datetime, timezone

    now = datetime.now(timezone.utc).isoformat()  # noqa: UP017

    for row in rows:
        pdm_id, pat_db_id, inv_id, patient_str = row
        sample_name = f"{patient_str}_sample_{pdm_id}"
        # Insert sample
        conn.execute(
            sa.text(
                "INSERT INTO samples (patient_id, sample_id, created_at) "
                "VALUES (:pat_id, :sid, :now)"
            ),
            {"pat_id": pat_db_id, "sid": sample_name, "now": now},
        )
        # Get the sample id
        sample_row = conn.execute(
            sa.text("SELECT id FROM samples WHERE sample_id = :sid"),
            {"sid": sample_name},
        ).fetchone()
        if sample_row:
            conn.execute(
                sa.text("UPDATE inventory_items SET sample_id = :sample_id WHERE id = :inv_id"),
                {"sample_id": sample_row[0], "inv_id": inv_id},
            )

    # ---- 4. Drop patient_data_map ----
    op.drop_table("patient_data_map")


def downgrade() -> None:
    # Recreate patient_data_map
    op.create_table(
        "patient_data_map",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "patient_id",
            sa.Integer(),
            sa.ForeignKey("patients.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "inventory_id",
            sa.Integer(),
            sa.ForeignKey("inventory_items.id", ondelete="CASCADE"),
            nullable=True,
        ),
        sa.Column("role", sa.String(), default="source"),
        sa.Column("notes", sa.Text(), nullable=True),
    )
    op.create_index("ix_patient_data_map_inventory_id", "patient_data_map", ["inventory_id"])

    # Migrate samples back to patient_data_map
    conn = op.get_bind()
    rows = conn.execute(
        sa.text(
            "SELECT s.patient_id, ii.id AS inventory_id "
            "FROM inventory_items ii "
            "JOIN samples s ON ii.sample_id = s.id "
            "WHERE s.patient_id IS NOT NULL"
        )
    ).fetchall()
    for row in rows:
        conn.execute(
            sa.text(
                "INSERT INTO patient_data_map (patient_id, inventory_id, role) "
                "VALUES (:pat_id, :inv_id, 'source')"
            ),
            {"pat_id": row[0], "inv_id": row[1]},
        )

    # Remove sample_id from inventory_items
    op.drop_index("ix_inventory_items_sample_id", "inventory_items")
    op.drop_column("inventory_items", "sample_id")

    # Drop samples table
    op.drop_index("ix_samples_patient_id", "samples")
    op.drop_table("samples")
