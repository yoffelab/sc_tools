"""Drop platform and domain columns from projects table.

Revision ID: 0015
Revises: 0014
Create Date: 2026-03-16

Platform and domain are properties of data sources, not projects.
Projects can span multiple platforms/domains via their linked datasets.
"""

from __future__ import annotations

revision = "0015"
down_revision = "0014"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    with op.batch_alter_table("projects") as batch_op:
        batch_op.drop_column("platform")
        batch_op.drop_column("domain")


def downgrade() -> None:
    with op.batch_alter_table("projects") as batch_op:
        batch_op.add_column(sa.Column("platform", sa.String(), nullable=True))
        batch_op.add_column(sa.Column("domain", sa.String(), nullable=True))
