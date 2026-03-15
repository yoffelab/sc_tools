"""Add modality column to bio_data table.

Revision ID: 0006
Revises: 0005
Create Date: 2026-03-11

Changes
-------
* Add ``modality`` column to ``bio_data`` -- denormalized from PlatformSpec
  for efficient queries without joining external Python registry.
* Backfill existing rows from PlatformSpec.modality where platform is known.
"""

from __future__ import annotations

revision = "0006"
down_revision = "0005"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    op.add_column("bio_data", sa.Column("modality", sa.String(), nullable=True))

    # Backfill from PlatformSpec registry
    try:
        from sc_tools.biodata import KNOWN_PLATFORMS

        conn = op.get_bind()
        for slug, spec in KNOWN_PLATFORMS.items():
            if spec.modality:
                conn.execute(
                    sa.text("UPDATE bio_data SET modality = :modality WHERE platform = :platform"),
                    {"modality": spec.modality, "platform": slug},
                )
    except ImportError:
        pass  # sc_tools not installed; skip backfill


def downgrade() -> None:
    op.drop_column("bio_data", "modality")
