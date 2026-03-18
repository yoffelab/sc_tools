"""Major schema redesign: collapse to 5 tables.

Revision ID: 0008
Revises: 0007
Create Date: 2026-03-13

Changes
-------
* Create ``patients`` table (replaces subjects + samples).
* Create ``data`` table (replaces bio_data + datasets + absorbs phase tracking).
* Migrate data from old tables into new ones.
* Drop 12 old tables: bio_images, rnaseq_data, spatial_seq_data,
  epigenomics_data, genome_seq_data, subject_project_links, samples,
  subjects, project_phases, datasets, bio_data, slurm_jobs.
* Simplify ``projects``: drop data_type, imaging_modality, project_type,
  visibility, phases_complete columns.
"""

from __future__ import annotations

revision = "0008"
down_revision = "0007"
branch_labels = None
depends_on = None

import sqlalchemy as sa
from alembic import op


def upgrade() -> None:
    # ---- 1. Create new tables ----

    op.create_table(
        "patients",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("patient_id", sa.String(), unique=True, nullable=False),
        sa.Column("metadata", sa.dialects.postgresql.JSONB(), server_default="{}"),
        sa.Column("created_at", sa.String()),
    )

    op.create_table(
        "data",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column(
            "project_id",
            sa.Integer(),
            sa.ForeignKey("projects.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column(
            "patient_id",
            sa.Integer(),
            sa.ForeignKey("patients.id", ondelete="SET NULL"),
            nullable=True,
        ),
        sa.Column("phase", sa.String()),
        sa.Column("status", sa.String(), server_default="pending"),
        sa.Column("uri", sa.Text(), nullable=False),
        sa.Column("format", sa.String()),
        sa.Column("platform", sa.String()),
        sa.Column("category", sa.String()),
        sa.Column("file_role", sa.String(), server_default="primary"),
        sa.Column("n_obs", sa.Integer()),
        sa.Column("n_vars", sa.Integer()),
        sa.Column("size_mb", sa.Float()),
        sa.Column("metadata", sa.dialects.postgresql.JSONB(), server_default="{}"),
        sa.Column("created_at", sa.String()),
        sa.Column("updated_at", sa.String()),
    )

    # ---- 2. Migrate data ----

    # 2a. subjects -> patients (merge sample info into metadata JSONB)
    op.execute(
        """
        INSERT INTO patients (id, patient_id, metadata, created_at)
        SELECT
            s.id,
            s.subject_id,
            jsonb_build_object(
                'organism', s.organism,
                'sex', s.sex,
                'age_at_collection', s.age_at_collection,
                'diagnosis', s.diagnosis,
                'diagnosis_code', s.diagnosis_code,
                'disease_stage', s.disease_stage,
                'treatment_status', s.treatment_status,
                'tissue_of_origin', s.tissue_of_origin,
                'cause_of_death', s.cause_of_death,
                'survival_days', s.survival_days,
                'vital_status', s.vital_status,
                'clinical_metadata', CASE
                    WHEN s.clinical_metadata_json IS NOT NULL
                    THEN s.clinical_metadata_json::jsonb
                    ELSE '{}'::jsonb
                END,
                'samples', COALESCE(
                    (SELECT jsonb_agg(jsonb_build_object(
                        'sample_id', sm.sample_id,
                        'tissue', sm.tissue,
                        'tissue_region', sm.tissue_region,
                        'fixation_method', sm.fixation_method,
                        'sample_type', sm.sample_type,
                        'batch', sm.batch,
                        'notes', sm.notes
                    ))
                    FROM samples sm WHERE sm.subject_id = s.id),
                    '[]'::jsonb
                )
            ),
            s.created_at
        FROM subjects s
        """
    )

    # Reset patients id sequence to max id
    op.execute(
        """
        SELECT setval(
            pg_get_serial_sequence('patients', 'id'),
            COALESCE((SELECT MAX(id) FROM patients), 1)
        )
        """
    )

    # 2b. bio_data -> data (merge JTI child columns into metadata JSONB)
    op.execute(
        """
        INSERT INTO data (
            id, project_id, patient_id, phase, status, uri, format,
            platform, category, file_role, n_obs, n_vars, size_mb,
            metadata, created_at, updated_at
        )
        SELECT
            bd.id,
            bd.project_id,
            -- Look up patient_id: bio_data.sample_id -> samples.subject_id -> patients.id
            (SELECT p.id FROM patients p
             JOIN subjects sub ON sub.subject_id = p.patient_id
             JOIN samples smp ON smp.subject_id = sub.id
             WHERE smp.id = bd.sample_id
             LIMIT 1),
            bd.phase,
            bd.status,
            bd.uri,
            bd.format,
            bd.platform,
            bd.category,
            bd.file_role,
            bd.n_obs,
            bd.n_vars,
            bd.size_mb,
            jsonb_build_object(
                'type', bd.type,
                'subcategory', bd.subcategory,
                'modality', bd.modality,
                'measurement', bd.measurement,
                'resolution', bd.resolution,
                'spatial', bd.spatial,
                'validated', bd.validated,
                'md5', bd.md5,
                'legacy_dataset_id', bd.legacy_dataset_id,
                'bio_image', (SELECT jsonb_build_object(
                    'image_type', bi.image_type,
                    'n_channels', bi.n_channels,
                    'pixel_size_um', bi.pixel_size_um,
                    'width_px', bi.width_px,
                    'height_px', bi.height_px,
                    'staining_protocol', bi.staining_protocol,
                    'channel_names_json', bi.channel_names_json
                ) FROM bio_images bi WHERE bi.id = bd.id),
                'spatial_seq', (SELECT jsonb_build_object(
                    'spatial_resolution', ss.spatial_resolution,
                    'panel_size', ss.panel_size,
                    'bin_size_um', ss.bin_size_um,
                    'fov_count', ss.fov_count,
                    'coordinate_system', ss.coordinate_system,
                    'tissue_area_mm2', ss.tissue_area_mm2,
                    'sequencing_platform', ss.sequencing_platform
                ) FROM spatial_seq_data ss WHERE ss.id = bd.id),
                'rnaseq', (SELECT jsonb_build_object(
                    'library_type', rn.library_type,
                    'sequencing_platform', rn.sequencing_platform,
                    'chemistry', rn.chemistry,
                    'read_length', rn.read_length,
                    'n_reads', rn.n_reads,
                    'reference_genome', rn.reference_genome,
                    'gene_annotation', rn.gene_annotation
                ) FROM rnaseq_data rn WHERE rn.id = bd.id),
                'epigenomics', (SELECT jsonb_build_object(
                    'assay_type', ep.assay_type,
                    'target_protein', ep.target_protein,
                    'n_peaks', ep.n_peaks,
                    'n_fragments', ep.n_fragments,
                    'genome_coverage', ep.genome_coverage
                ) FROM epigenomics_data ep WHERE ep.id = bd.id),
                'genome_seq', (SELECT jsonb_build_object(
                    'sequencing_type', gs.sequencing_type,
                    'panel_name', gs.panel_name,
                    'coverage_mean', gs.coverage_mean,
                    'reference_genome', gs.reference_genome,
                    'sequencing_platform', gs.sequencing_platform
                ) FROM genome_seq_data gs WHERE gs.id = bd.id)
            ),
            bd.created_at,
            bd.updated_at
        FROM bio_data bd
        """
    )

    # Reset data id sequence to max id
    op.execute(
        """
        SELECT setval(
            pg_get_serial_sequence('data', 'id'),
            COALESCE((SELECT MAX(id) FROM data), 1)
        )
        """
    )

    # ---- 3. Drop old tables ----
    # Must drop in correct FK dependency order.
    # JTI children first (reference bio_data):
    op.drop_table("bio_images")
    op.drop_table("rnaseq_data")
    op.drop_table("spatial_seq_data")
    op.drop_table("epigenomics_data")
    op.drop_table("genome_seq_data")
    # project_phases references bio_data
    op.drop_table("project_phases")
    # datasets has FK to bio_data and bio_data has FK to datasets
    # Drop the bio_data -> datasets FK first, then drop both
    op.drop_constraint("bio_data_legacy_dataset_id_fkey", "bio_data", type_="foreignkey")
    op.drop_constraint("bio_data_sample_id_fkey", "bio_data", type_="foreignkey")
    op.drop_constraint("bio_data_subject_id_fkey", "bio_data", type_="foreignkey")
    op.drop_table("datasets")
    op.drop_table("bio_data")
    # Now safe to drop subject/sample tables
    op.drop_table("subject_project_links")
    op.drop_table("samples")
    op.drop_table("subjects")
    op.drop_table("slurm_jobs")

    # ---- 4. Simplify projects table ----

    op.drop_column("projects", "data_type")
    op.drop_column("projects", "imaging_modality")
    op.drop_column("projects", "project_type")
    op.drop_column("projects", "visibility")
    op.drop_column("projects", "phases_complete")


def downgrade() -> None:
    raise NotImplementedError(
        "Downgrade not supported for migration 0008 (major schema redesign). "
        "Restore from a database backup."
    )
