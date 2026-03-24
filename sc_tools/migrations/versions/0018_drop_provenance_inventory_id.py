"""Drop inventory_id from provenance table.

Revision ID: 0018
Revises: 0017
Create Date: 2026-03-17

Inventory items are ingest artifacts, not transformation outputs.
Provenance only applies to phases (transformations) and datasets (assembly).
"""

from __future__ import annotations

revision = "0018"
down_revision = "0017"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # Delete any provenance rows targeting inventory items
    conn = op.get_bind()
    conn.execute(sa.text("DELETE FROM provenance WHERE target_type = 'inventory'"))

    # Drop the index and column
    op.drop_index("ix_provenance_inventory", "provenance")
    op.drop_column("provenance", "inventory_id")


def downgrade() -> None:
    op.add_column(
        "provenance",
        sa.Column(
            "inventory_id",
            sa.Integer(),
            sa.ForeignKey("inventory_items.id", ondelete="SET NULL"),
            nullable=True,
        ),
    )
    op.create_index("ix_provenance_inventory", "provenance", ["inventory_id"])
