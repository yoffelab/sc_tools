"""Add data_sources and project_data_sources tables.

Revision ID: 0004
Revises: 0003
Create Date: 2026-03-06

Changes
-------
* New ``data_sources`` table — raw data catalog (HPC directories, public
  datasets, external repositories).  Distinct from ``projects`` (lab work)
  and ``datasets`` (processed checkpoints).
* New ``project_data_sources`` join table — many-to-many between projects
  and data sources.

source_type vocabulary
----------------------
    hpc_lab            -- data produced by our lab on HPC scratch
    hpc_collaborator   -- collaborator data on shared HPC storage
    public_10x         -- 10x Genomics public datasets portal
    public_geo         -- NCBI GEO accession
    public_zenodo      -- Zenodo DOI
    public_portal      -- other public portal (HuBMAP, HTAN, HCA, etc.)

status vocabulary
-----------------
    available          -- data is accessible at the given URI
    restricted         -- requires special access / VPN / permission
    pending_download   -- URL known but not yet downloaded
    archived           -- moved to long-term storage
"""

from __future__ import annotations

revision = "0004"
down_revision = "0003"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    op.create_table(
        "data_sources",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("name", sa.String(), unique=True, nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("uri", sa.Text(), nullable=False),
        sa.Column("platform", sa.String(), nullable=True),
        sa.Column("domain", sa.String(), nullable=True),
        sa.Column("imaging_modality", sa.String(), nullable=True),
        sa.Column("source_type", sa.String(), nullable=True),
        sa.Column("organism", sa.String(), nullable=True),
        sa.Column("tissue", sa.String(), nullable=True),
        sa.Column("disease", sa.String(), nullable=True),
        sa.Column("n_samples", sa.Integer(), nullable=True),
        sa.Column("n_cells", sa.Integer(), nullable=True),
        sa.Column("publication", sa.String(), nullable=True),
        sa.Column("access_notes", sa.Text(), nullable=True),
        sa.Column("status", sa.String(), server_default="available", nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
    )

    op.create_table(
        "project_data_sources",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "data_source_id",
            sa.Integer(),
            sa.ForeignKey("data_sources.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("role", sa.String(), server_default="input", nullable=True),
        sa.Column("notes", sa.Text(), nullable=True),
    )


def downgrade() -> None:
    op.drop_table("project_data_sources")
    op.drop_table("data_sources")
