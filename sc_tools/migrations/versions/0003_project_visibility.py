"""Project visibility and access control.

Revision ID: 0003
Revises: 0002
Create Date: 2026-03-06

Changes
-------
* ``projects`` table: add ``project_type`` and ``visibility`` columns.

project_type vocabulary
------------------------
    internal    -- lab-led project (default)
    external    -- collaboration or contract project

visibility vocabulary
----------------------
    private     -- restricted access; data not publicly shareable (default)
    public      -- shareable / published; data may be deposited or open-access
"""

from __future__ import annotations

revision = "0003"
down_revision = "0002"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    with op.batch_alter_table("projects") as batch_op:
        batch_op.add_column(
            sa.Column("project_type", sa.String(), nullable=True, server_default="internal")
        )
        batch_op.add_column(
            sa.Column("visibility", sa.String(), nullable=True, server_default="private")
        )


def downgrade() -> None:
    with op.batch_alter_table("projects") as batch_op:
        batch_op.drop_column("visibility")
        batch_op.drop_column("project_type")
