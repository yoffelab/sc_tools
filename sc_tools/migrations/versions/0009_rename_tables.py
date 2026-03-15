"""Rename data → data_phase, data_sources → data_inventory, project_data_sources → project_data_inventory.

Revision ID: 0009
Revises: 0008
Create Date: 2026-03-13
"""

from __future__ import annotations

revision = "0009"
down_revision = "0008"
branch_labels = None
depends_on = None

from alembic import op


def upgrade() -> None:
    op.rename_table("data", "data_phase")
    op.rename_table("data_sources", "data_inventory")
    op.rename_table("project_data_sources", "project_data_inventory")


def downgrade() -> None:
    op.rename_table("project_data_inventory", "project_data_sources")
    op.rename_table("data_inventory", "data_sources")
    op.rename_table("data_phase", "data")
