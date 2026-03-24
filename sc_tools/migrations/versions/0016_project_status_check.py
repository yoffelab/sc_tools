"""Add CHECK constraint on projects.status.

Revision ID: 0016
Revises: 0015
Create Date: 2026-03-16

Valid values: active, paused, completed, archived.
"""

from __future__ import annotations

revision = "0016"
down_revision = "0015"
branch_labels = None
depends_on = None

from alembic import op


def upgrade() -> None:
    op.create_check_constraint(
        "ck_projects_status",
        "projects",
        "status IN ('active', 'paused', 'completed', 'archived')",
    )


def downgrade() -> None:
    op.drop_constraint("ck_projects_status", "projects", type_="check")
