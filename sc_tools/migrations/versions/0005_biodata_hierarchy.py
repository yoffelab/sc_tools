"""Add BioData type hierarchy, subjects, and samples tables.

Revision ID: 0005
Revises: 0004
Create Date: 2026-03-11

Changes
-------
* New ``subjects`` table -- cross-project de-identified clinical metadata.
* New ``samples`` table -- links subjects to physical specimens within projects.
* New ``subject_project_links`` join table -- subjects can appear across projects.
* New ``bio_data`` table -- JTI base for typed biological data objects.
* New ``bio_images`` table -- image-specific columns (JTI child of bio_data).
* New ``rnaseq_data`` table -- RNA-seq-specific columns.
* New ``spatial_seq_data`` table -- spatial sequencing-specific columns.
* New ``epigenomics_data`` table -- epigenomics-specific columns.
* New ``genome_seq_data`` table -- genome sequencing-specific columns.
* Add ``bio_data_id`` FK on ``datasets`` for forward-linking during migration.
* The ``datasets`` table is NOT dropped -- backward compatibility preserved.
"""

from __future__ import annotations

revision = "0005"
down_revision = "0004"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # -- subjects --
    op.create_table(
        "subjects",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("subject_id", sa.String(), unique=True, nullable=False),
        sa.Column("organism", sa.String(), server_default="human"),
        sa.Column("sex", sa.String(), nullable=True),
        sa.Column("age_at_collection", sa.Float(), nullable=True),
        sa.Column("diagnosis", sa.String(), nullable=True),
        sa.Column("diagnosis_code", sa.String(), nullable=True),
        sa.Column("disease_stage", sa.String(), nullable=True),
        sa.Column("treatment_status", sa.String(), nullable=True),
        sa.Column("tissue_of_origin", sa.String(), nullable=True),
        sa.Column("cause_of_death", sa.String(), nullable=True),
        sa.Column("survival_days", sa.Float(), nullable=True),
        sa.Column("vital_status", sa.String(), nullable=True),
        sa.Column("clinical_metadata_json", sa.Text(), nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
    )

    # -- samples --
    op.create_table(
        "samples",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("sample_id", sa.String(), nullable=False),
        sa.Column(
            "subject_id",
            sa.Integer(),
            sa.ForeignKey("subjects.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=True,
        ),
        sa.Column("tissue", sa.String(), nullable=True),
        sa.Column("tissue_region", sa.String(), nullable=True),
        sa.Column("collection_date", sa.String(), nullable=True),
        sa.Column("fixation_method", sa.String(), nullable=True),
        sa.Column("sample_type", sa.String(), nullable=True),
        sa.Column("batch", sa.String(), nullable=True),
        sa.Column("notes", sa.Text(), nullable=True),
        sa.Column("created_at", sa.String(), nullable=True),
    )

    # -- subject_project_links --
    op.create_table(
        "subject_project_links",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "subject_id",
            sa.Integer(),
            sa.ForeignKey("subjects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("role", sa.String(), server_default="enrolled"),
        sa.Column("notes", sa.Text(), nullable=True),
        sa.UniqueConstraint("subject_id", "project_id", name="uq_subject_project"),
    )

    # -- bio_data (JTI base) --
    op.create_table(
        "bio_data",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("type", sa.String(), nullable=False),  # JTI discriminator
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "subject_id",
            sa.Integer(),
            sa.ForeignKey("subjects.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column(
            "sample_id",
            sa.Integer(),
            sa.ForeignKey("samples.id", ondelete="SET NULL"),
            nullable=True,
        ),
        # classification
        sa.Column("category", sa.String(), nullable=True),
        sa.Column("subcategory", sa.String(), nullable=True),
        sa.Column("platform", sa.String(), nullable=True),
        sa.Column("measurement", sa.String(), nullable=True),
        sa.Column("resolution", sa.String(), nullable=True),
        sa.Column("spatial", sa.Boolean(), nullable=True),
        # data tracking
        sa.Column("uri", sa.Text(), nullable=False),
        sa.Column("format", sa.String(), nullable=True),
        sa.Column("status", sa.String(), server_default="pending"),
        sa.Column("file_role", sa.String(), server_default="primary"),
        sa.Column("validated", sa.Boolean(), server_default="0"),
        sa.Column("n_obs", sa.Integer(), nullable=True),
        sa.Column("n_vars", sa.Integer(), nullable=True),
        sa.Column("phase", sa.String(), nullable=True),
        sa.Column("size_mb", sa.Float(), nullable=True),
        sa.Column("md5", sa.String(), nullable=True),
        # provenance
        sa.Column(
            "legacy_dataset_id",
            sa.Integer(),
            sa.ForeignKey("datasets.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column("created_at", sa.String(), nullable=True),
        sa.Column("updated_at", sa.String(), nullable=True),
    )

    # -- bio_images (JTI child) --
    op.create_table(
        "bio_images",
        sa.Column(
            "id",
            sa.Integer(),
            sa.ForeignKey("bio_data.id", ondelete="CASCADE"),
            primary_key=True,
        ),
        sa.Column("image_type", sa.String(), nullable=True),
        sa.Column("n_channels", sa.Integer(), nullable=True),
        sa.Column("pixel_size_um", sa.Float(), nullable=True),
        sa.Column("width_px", sa.Integer(), nullable=True),
        sa.Column("height_px", sa.Integer(), nullable=True),
        sa.Column("staining_protocol", sa.String(), nullable=True),
        sa.Column("channel_names_json", sa.Text(), nullable=True),
    )

    # -- rnaseq_data (JTI child) --
    op.create_table(
        "rnaseq_data",
        sa.Column(
            "id",
            sa.Integer(),
            sa.ForeignKey("bio_data.id", ondelete="CASCADE"),
            primary_key=True,
        ),
        sa.Column("library_type", sa.String(), nullable=True),
        sa.Column("sequencing_platform", sa.String(), nullable=True),
        sa.Column("chemistry", sa.String(), nullable=True),
        sa.Column("read_length", sa.Integer(), nullable=True),
        sa.Column("n_reads", sa.Integer(), nullable=True),
        sa.Column("reference_genome", sa.String(), nullable=True),
        sa.Column("gene_annotation", sa.String(), nullable=True),
    )

    # -- spatial_seq_data (JTI child) --
    op.create_table(
        "spatial_seq_data",
        sa.Column(
            "id",
            sa.Integer(),
            sa.ForeignKey("bio_data.id", ondelete="CASCADE"),
            primary_key=True,
        ),
        sa.Column("spatial_resolution", sa.String(), nullable=True),
        sa.Column("panel_size", sa.Integer(), nullable=True),
        sa.Column("bin_size_um", sa.Float(), nullable=True),
        sa.Column("fov_count", sa.Integer(), nullable=True),
        sa.Column("coordinate_system", sa.String(), nullable=True),
        sa.Column("tissue_area_mm2", sa.Float(), nullable=True),
        sa.Column("sequencing_platform", sa.String(), nullable=True),
    )

    # -- epigenomics_data (JTI child) --
    op.create_table(
        "epigenomics_data",
        sa.Column(
            "id",
            sa.Integer(),
            sa.ForeignKey("bio_data.id", ondelete="CASCADE"),
            primary_key=True,
        ),
        sa.Column("assay_type", sa.String(), nullable=True),
        sa.Column("target_protein", sa.String(), nullable=True),
        sa.Column("n_peaks", sa.Integer(), nullable=True),
        sa.Column("n_fragments", sa.Integer(), nullable=True),
        sa.Column("genome_coverage", sa.Float(), nullable=True),
    )

    # -- genome_seq_data (JTI child) --
    op.create_table(
        "genome_seq_data",
        sa.Column(
            "id",
            sa.Integer(),
            sa.ForeignKey("bio_data.id", ondelete="CASCADE"),
            primary_key=True,
        ),
        sa.Column("sequencing_type", sa.String(), nullable=True),
        sa.Column("panel_name", sa.String(), nullable=True),
        sa.Column("coverage_mean", sa.Float(), nullable=True),
        sa.Column("reference_genome", sa.String(), nullable=True),
        sa.Column("sequencing_platform", sa.String(), nullable=True),
    )

    # -- Add bio_data_id FK on datasets for forward-linking --
    op.add_column(
        "datasets",
        sa.Column(
            "bio_data_id",
            sa.Integer(),
            sa.ForeignKey("bio_data.id", ondelete="SET NULL"),
            nullable=True,
        ),
    )


def downgrade() -> None:
    op.drop_column("datasets", "bio_data_id")
    op.drop_table("genome_seq_data")
    op.drop_table("epigenomics_data")
    op.drop_table("spatial_seq_data")
    op.drop_table("rnaseq_data")
    op.drop_table("bio_images")
    op.drop_table("bio_data")
    op.drop_table("subject_project_links")
    op.drop_table("samples")
    op.drop_table("subjects")
