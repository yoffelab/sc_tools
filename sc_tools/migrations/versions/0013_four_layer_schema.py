"""Add four-layer schema: data_sources, inventory_items, datasets, project_phases, provenance.

Revision ID: 0013
Revises: 0012
Create Date: 2026-03-16

Additive only -- old and new tables coexist. No drops, no renames of existing tables.

New tables
----------
* ``data_sources`` -- raw data sources (copy of data_inventory + metadata column)
* ``inventory_items`` -- clean ingested data (h5ad files)
* ``datasets`` -- versioned assemblies (MuData or AnnData)
* ``dataset_members`` -- composition of a dataset
* ``project_datasets`` -- project-to-dataset links (replaces data_project_map)
* ``project_phases`` -- phase tracking per project+dataset
* ``provenance`` -- tool/version/params tracking

Data migration
--------------
* Copy data_inventory -> data_sources
* Migrate data_processing_phase rows with real URIs -> inventory_items
* Migrate phase markers -> project_phases
* Add inventory_id column to patient_data_map
* Copy data_project_map links -> project_datasets (with placeholder datasets)
"""

from __future__ import annotations

revision = "0013"
down_revision = "0012"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # ---- 1. Create new tables ----

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
        sa.Column("status", sa.String(), server_default="discovered", nullable=True),
        sa.Column("metadata", sa.JSON(), server_default="{}", nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
    )

    op.create_table(
        "inventory_items",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "data_source_id",
            sa.Integer(),
            sa.ForeignKey("data_sources.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column("name", sa.String(), unique=True, nullable=False),
        sa.Column("uri", sa.Text(), nullable=False),
        sa.Column("modality", sa.String(), nullable=False),
        sa.Column("platform", sa.String(), nullable=True),
        sa.Column("format", sa.String(), server_default="h5ad", nullable=True),
        sa.Column("n_obs", sa.Integer(), nullable=True),
        sa.Column("n_vars", sa.Integer(), nullable=True),
        sa.Column("size_mb", sa.Float(), nullable=True),
        sa.Column("organism", sa.String(), nullable=True),
        sa.Column("tissue", sa.String(), nullable=True),
        sa.Column("metadata", sa.JSON(), server_default="{}", nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
        sa.Column("updated_at", sa.String(), nullable=True),
    )

    op.create_table(
        "datasets",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("name", sa.String(), nullable=False),
        sa.Column("version", sa.Integer(), nullable=False, server_default="1"),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("uri", sa.Text(), nullable=True),
        sa.Column("format", sa.String(), server_default="mudata", nullable=True),
        sa.Column("n_obs", sa.Integer(), nullable=True),
        sa.Column("size_mb", sa.Float(), nullable=True),
        sa.Column("is_current", sa.Boolean(), server_default=sa.text("true"), nullable=True),
        sa.Column("metadata", sa.JSON(), server_default="{}", nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
        sa.Column("updated_at", sa.String(), nullable=True),
        sa.UniqueConstraint("name", "version", name="uq_datasets_name_version"),
    )

    op.create_table(
        "dataset_members",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "dataset_id",
            sa.Integer(),
            sa.ForeignKey("datasets.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "inventory_id",
            sa.Integer(),
            sa.ForeignKey("inventory_items.id", ondelete="RESTRICT"),
            nullable=False,
        ),
        sa.Column("modality_key", sa.String(), nullable=False),
        sa.UniqueConstraint(
            "dataset_id", "modality_key", name="uq_dataset_members_dataset_modality"
        ),
    )

    op.create_table(
        "project_datasets",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "dataset_id",
            sa.Integer(),
            sa.ForeignKey("datasets.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("role", sa.String(), server_default="primary", nullable=True),
        sa.Column("notes", sa.Text(), nullable=True),
        sa.UniqueConstraint("project_id", "dataset_id", name="uq_project_datasets_project_dataset"),
    )

    op.create_table(
        "project_phases",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "dataset_id",
            sa.Integer(),
            sa.ForeignKey("datasets.id", ondelete="RESTRICT"),
            nullable=False,
        ),
        sa.Column("phase_group", sa.String(), nullable=False),
        sa.Column("subphase", sa.String(), nullable=False),
        sa.Column("status", sa.String(), server_default="pending", nullable=True),
        sa.Column("uri", sa.Text(), nullable=True),
        sa.Column("n_obs", sa.Integer(), nullable=True),
        sa.Column("n_vars", sa.Integer(), nullable=True),
        sa.Column("metadata", sa.JSON(), server_default="{}", nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
        sa.Column("updated_at", sa.String(), nullable=True),
        sa.UniqueConstraint(
            "project_id",
            "dataset_id",
            "phase_group",
            "subphase",
            name="uq_project_phases_project_dataset_group_subphase",
        ),
        sa.CheckConstraint(
            "phase_group IN ('data_processing', 'discovery')",
            name="ck_project_phases_phase_group",
        ),
    )

    op.create_table(
        "provenance",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("target_type", sa.String(), nullable=False),
        sa.Column(
            "inventory_id",
            sa.Integer(),
            sa.ForeignKey("inventory_items.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column(
            "phase_id",
            sa.Integer(),
            sa.ForeignKey("project_phases.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column(
            "dataset_id",
            sa.Integer(),
            sa.ForeignKey("datasets.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column("tool", sa.String(), nullable=False),
        sa.Column("tool_version", sa.String(), nullable=True),
        sa.Column("reference_genome", sa.String(), nullable=True),
        sa.Column("reference_dataset", sa.Text(), nullable=True),
        sa.Column("signature_source", sa.Text(), nullable=True),
        sa.Column("n_input_obs", sa.Integer(), nullable=True),
        sa.Column("n_output_obs", sa.Integer(), nullable=True),
        sa.Column("params", sa.JSON(), server_default="{}", nullable=True),
        sa.Column("environment", sa.JSON(), server_default="{}", nullable=True),
        sa.Column("script_uri", sa.Text(), nullable=True),
        sa.Column("agent", sa.String(), nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
        sa.CheckConstraint(
            "(target_type = 'inventory' AND inventory_id IS NOT NULL AND phase_id IS NULL AND dataset_id IS NULL) OR "
            "(target_type = 'phase' AND inventory_id IS NULL AND phase_id IS NOT NULL AND dataset_id IS NULL) OR "
            "(target_type = 'dataset' AND inventory_id IS NULL AND phase_id IS NULL AND dataset_id IS NOT NULL)",
            name="ck_provenance_target_type",
        ),
    )

    # ---- 2. Create indexes ----

    op.create_index("ix_inventory_items_data_source_id", "inventory_items", ["data_source_id"])
    op.create_index("ix_dataset_members_dataset_id", "dataset_members", ["dataset_id"])
    op.create_index("ix_project_datasets_project_id", "project_datasets", ["project_id"])
    op.create_index("ix_project_phases_project_id", "project_phases", ["project_id"])
    op.create_index(
        "ix_project_phases_project_dataset", "project_phases", ["project_id", "dataset_id"]
    )
    op.create_index("ix_provenance_inventory", "provenance", ["inventory_id"])
    op.create_index("ix_provenance_phase", "provenance", ["phase_id"])
    op.create_index("ix_provenance_dataset", "provenance", ["dataset_id"])

    # ---- 3. Copy data_inventory -> data_sources ----

    op.execute(
        """
        INSERT INTO data_sources (
            id, name, description, uri, platform, domain, imaging_modality,
            source_type, organism, tissue, disease, n_samples, n_cells,
            publication, access_notes, status, created_at
        )
        SELECT
            id, name, description, uri, platform, domain, imaging_modality,
            source_type, organism, tissue, disease, n_samples, n_cells,
            publication, access_notes, status, created_at
        FROM data_inventory
        """
    )

    # ---- 4. Migrate data_processing_phase rows with real URIs -> inventory_items ----
    # Real file rows: uri does NOT start with 'phase://' and file_role != 'phase_marker'

    op.execute(
        """
        INSERT INTO inventory_items (
            id, name, uri, modality, platform, format,
            n_obs, n_vars, size_mb, metadata, created_at, updated_at
        )
        SELECT
            id,
            COALESCE(phase, 'unknown') || '_' || CAST(id AS VARCHAR),
            uri,
            COALESCE(category, 'unknown'),
            platform,
            format,
            n_obs,
            n_vars,
            size_mb,
            metadata,
            created_at,
            updated_at
        FROM data_processing_phase
        WHERE uri NOT LIKE 'phase://%%'
          AND COALESCE(file_role, 'primary') != 'phase_marker'
        """
    )

    # ---- 5. Migrate phase markers -> project_phases ----
    # Create placeholder datasets for each project (needed for project_phases FK)

    op.execute(
        """
        INSERT INTO datasets (name, version, description, format, is_current, created_at)
        SELECT
            'legacy_' || p.name,
            1,
            'Auto-created placeholder for migrated phase markers from project ' || p.name,
            'h5ad',
            TRUE,
            p.created_at
        FROM projects p
        WHERE EXISTS (
            SELECT 1 FROM data_processing_phase dpp
            WHERE dpp.project_id = p.id
              AND (dpp.file_role = 'phase_marker' OR dpp.uri LIKE 'phase://%%')
        )
        """
    )

    # Now insert phase markers into project_phases
    op.execute(
        """
        INSERT INTO project_phases (
            project_id, dataset_id, phase_group, subphase, status,
            uri, n_obs, n_vars, metadata, created_at, updated_at
        )
        SELECT
            dpp.project_id,
            d.id,
            'data_processing',
            dpp.phase,
            dpp.status,
            dpp.uri,
            dpp.n_obs,
            dpp.n_vars,
            dpp.metadata,
            dpp.created_at,
            dpp.updated_at
        FROM data_processing_phase dpp
        JOIN projects p ON p.id = dpp.project_id
        JOIN datasets d ON d.name = 'legacy_' || p.name
        WHERE dpp.file_role = 'phase_marker' OR dpp.uri LIKE 'phase://%%'
        """
    )

    # ---- 6. Add inventory_id column to patient_data_map ----

    with op.batch_alter_table("patient_data_map") as batch_op:
        batch_op.add_column(sa.Column("inventory_id", sa.Integer(), nullable=True))

    # Populate inventory_id from data_id where the data row was migrated
    op.execute(
        """
        UPDATE patient_data_map
        SET inventory_id = data_id
        WHERE data_id IN (SELECT id FROM inventory_items)
        """
    )

    op.create_index("ix_patient_data_map_inventory_id", "patient_data_map", ["inventory_id"])

    # ---- 7. Copy data_project_map links -> project_datasets ----
    # Create placeholder datasets for projects that have data_project_map links
    # but no legacy dataset yet

    op.execute(
        """
        INSERT INTO datasets (name, version, description, format, is_current, created_at)
        SELECT
            'legacy_' || p.name,
            1,
            'Auto-created placeholder for project ' || p.name,
            'h5ad',
            TRUE,
            p.created_at
        FROM projects p
        JOIN data_project_map dpm ON dpm.project_id = p.id
        WHERE NOT EXISTS (
            SELECT 1 FROM datasets d WHERE d.name = 'legacy_' || p.name
        )
        GROUP BY p.id, p.name, p.created_at
        """
    )

    op.execute(
        """
        INSERT INTO project_datasets (project_id, dataset_id, role, notes)
        SELECT DISTINCT ON (dpm.project_id, d.id)
            dpm.project_id,
            d.id,
            dpm.role,
            dpm.notes
        FROM data_project_map dpm
        JOIN projects p ON p.id = dpm.project_id
        JOIN datasets d ON d.name = 'legacy_' || p.name
        ORDER BY dpm.project_id, d.id, dpm.id
        """
    )


def downgrade() -> None:
    # Drop indexes on patient_data_map
    op.drop_index("ix_patient_data_map_inventory_id", table_name="patient_data_map")

    # Remove inventory_id column from patient_data_map
    with op.batch_alter_table("patient_data_map") as batch_op:
        batch_op.drop_column("inventory_id")

    # Drop new tables in reverse dependency order
    op.drop_table("provenance")
    op.drop_table("project_phases")
    op.drop_table("project_datasets")
    op.drop_table("dataset_members")
    op.drop_table("datasets")
    op.drop_table("inventory_items")
    op.drop_table("data_sources")
