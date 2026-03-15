"""Rename data_phase → data_processing_phase.

Revision ID: 0011
Revises: 0010
Create Date: 2026-03-13
"""

from __future__ import annotations

revision = "0011"
down_revision = "0010"
branch_labels = None
depends_on = None

from alembic import op


def upgrade() -> None:
    op.rename_table("data_phase", "data_processing_phase")


def downgrade() -> None:
    op.rename_table("data_processing_phase", "data_phase")
