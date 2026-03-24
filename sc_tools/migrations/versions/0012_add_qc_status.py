"""Add qc_status column to data_processing_phase.

Revision ID: 0012
Revises: 0011
Create Date: 2026-03-14

Values: pending | pass | fail | na
"""

from __future__ import annotations

revision = "0012"
down_revision = "0011"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    op.add_column(
        "data_processing_phase",
        sa.Column("qc_status", sa.String(), server_default="pending", nullable=True),
    )
    # Backfill existing data with status=ready as qc pass
    op.execute("UPDATE data_processing_phase SET qc_status = 'pass' WHERE status = 'ready'")


def downgrade() -> None:
    op.drop_column("data_processing_phase", "qc_status")
