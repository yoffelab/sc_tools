"""Data and job registry for sc_tools projects.

Tracks datasets (AnnData checkpoints), SLURM jobs, and agent tasks across
all active projects.  Uses **SQLite** by default (zero-config) with optional
**PostgreSQL** support via the ``SC_TOOLS_REGISTRY_URL`` environment variable.

Quick start
-----------
::

    from sc_tools.registry import Registry

    reg = Registry()
    reg.add_project("ggo_visium", platform="visium", data_type="visium",
                    domain="spatial_transcriptomics", imaging_modality="sequencing_based")
    ds_id = reg.register_dataset(
        "ggo_visium", phase="qc_filter",
        uri="sftp://brb//athena/.../adata.raw.h5ad",
        fmt="h5ad",
        file_role="primary",
    )
    reg.upsert_phase("ggo_visium", "qc_filter", status="complete",
                     primary_dataset_id=ds_id, n_obs=50000, n_samples=8)

CLI
---
::

    python -m sc_tools registry status

Environment
-----------
``SC_TOOLS_REGISTRY_URL``
    Optional SQLAlchemy-compatible DB URL.  Defaults to
    ``sqlite:///~/.sc_tools/registry.db``.

Technology taxonomy
-------------------
domain:
    ``spatial_transcriptomics`` | ``spatial_proteomics`` | ``imaging`` |
    ``single_cell`` | ``bulk``

imaging_modality:
    ``brightfield`` | ``fluorescence`` | ``multiplexed_fluorescence`` |
    ``probe_based`` | ``mass_spec_imaging`` | ``sequencing_based``

Phase slugs (new names → old codes)
------------------------------------
    ingest_raw       p0a
    ingest_load      p0b
    qc_filter        p1
    metadata_attach  p2
    preprocess       p3
    demographics     p3.5
    scoring          p3.5b
    celltype_manual  p4
    biology          p5
    meta_analysis    p6/p7
"""

from __future__ import annotations

import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_DEFAULT_DB_PATH = Path.home() / ".sc_tools" / "registry.db"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_db_url() -> str:
    """Return the database URL from env or default SQLite path."""
    url = os.environ.get("SC_TOOLS_REGISTRY_URL")
    if url:
        return url
    db_path = _DEFAULT_DB_PATH
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return f"sqlite:///{db_path}"


def _require_sqlalchemy() -> Any:
    try:
        import sqlalchemy

        return sqlalchemy
    except ImportError as e:
        raise ImportError(
            "SQLAlchemy is required for the registry. Install with: pip install sc-tools[registry]"
        ) from e


def _utcnow() -> str:
    return datetime.now(tz=timezone.utc).isoformat()  # noqa: UP017


# ---------------------------------------------------------------------------
# ORM models — built lazily to avoid top-level SQLAlchemy import errors
# ---------------------------------------------------------------------------


def _build_models(Base: Any) -> tuple:
    """Return ORM model classes.

    Returns a 16-tuple:
        Project, Dataset, SlurmJob, AgentTask, ProjectPhase, DataSource,
        ProjectDataSource, Subject, Sample, SubjectProjectLink, BioData,
        BioImage, RNASeqData, SpatialSeqData, EpigenomicsData, GenomeSeqData
    """
    from sqlalchemy import (
        Boolean,
        Column,
        Float,
        ForeignKey,
        Integer,
        String,
        Text,
        UniqueConstraint,
    )
    from sqlalchemy.orm import relationship

    class Project(Base):
        __tablename__ = "projects"

        id = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String, unique=True, nullable=False)
        platform = Column(String)
        data_type = Column(String)
        # Technology taxonomy (added in migration 0002)
        domain = Column(
            String
        )  # spatial_transcriptomics | spatial_proteomics | imaging | single_cell | bulk
        imaging_modality = Column(
            String
        )  # brightfield | fluorescence | multiplexed_fluorescence | probe_based | mass_spec_imaging | sequencing_based
        # Visibility / access control (added in migration 0003)
        project_type = Column(String, default="internal")  # internal | external
        visibility = Column(String, default="private")  # private | public
        phases_complete = Column(
            Text, default="[]"
        )  # JSON array (legacy; use project_phases table)
        status = Column(String, default="active")
        created_at = Column(String, default=_utcnow)

        datasets = relationship("Dataset", back_populates="project", cascade="all, delete-orphan")
        jobs = relationship("SlurmJob", back_populates="project", cascade="all, delete-orphan")
        tasks = relationship("AgentTask", back_populates="project", cascade="all, delete-orphan")
        phases = relationship(
            "ProjectPhase", back_populates="project", cascade="all, delete-orphan"
        )

    class Dataset(Base):
        __tablename__ = "datasets"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(Integer, ForeignKey("projects.id"), nullable=False)
        sample_id = Column(String)
        phase = Column(String)  # slug: qc_filter, scoring, celltype_manual, etc.
        uri = Column(Text, nullable=False)
        format = Column(String)  # h5ad, zarr, tiff, tsv
        size_mb = Column(Float)
        md5 = Column(String)
        status = Column(String, default="pending")  # pending, ready, archived, error
        # Extended metadata (added in migration 0002)
        file_role = Column(
            String, default="primary"
        )  # primary | supplementary | entry_point | spatialdata | image | metadata
        validated = Column(Boolean, default=False)
        n_obs = Column(Integer)
        n_vars = Column(Integer)
        # Forward-link to BioData (migration 0005)
        bio_data_id = Column(Integer, ForeignKey("bio_data.id", ondelete="SET NULL"), nullable=True)
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        project = relationship("Project", back_populates="datasets")

    class SlurmJob(Base):
        __tablename__ = "slurm_jobs"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(Integer, ForeignKey("projects.id"), nullable=False)
        sample_id = Column(String)
        phase = Column(String)
        cluster = Column(String)  # brb, cayuga
        slurm_job_id = Column(String)
        status = Column(String, default="submitted")  # submitted, running, completed, failed
        submitted_at = Column(String)
        finished_at = Column(String)
        log_uri = Column(Text)
        error_msg = Column(Text)

        project = relationship("Project", back_populates="jobs")

    class AgentTask(Base):
        __tablename__ = "agent_tasks"

        id = Column(Integer, primary_key=True, autoincrement=True)
        task_type = Column(String)
        project_id = Column(Integer, ForeignKey("projects.id"), nullable=False)
        status = Column(String, default="queued")  # queued, running, completed, failed
        inputs_json = Column(Text)
        outputs_json = Column(Text)
        error = Column(Text)
        started_at = Column(String)
        finished_at = Column(String)

        project = relationship("Project", back_populates="tasks")

    class ProjectPhase(Base):
        """Per-phase pipeline status row.

        Composite primary key ``(project_id, phase)``.  Upserted (not inserted)
        on each status update.  Use :meth:`Registry.upsert_phase` to write.
        """

        __tablename__ = "project_phases"

        project_id = Column(
            Integer, ForeignKey("projects.id", ondelete="CASCADE"), primary_key=True
        )
        phase = Column(String, primary_key=True)  # slug, e.g. "qc_filter"
        status = Column(
            String, default="not_started"
        )  # not_started | in_progress | complete | failed | skipped
        entry_phase = Column(Boolean, default=False)  # True = where the pipeline was loaded from
        primary_dataset_id = Column(
            Integer, ForeignKey("datasets.id"), nullable=True
        )  # FK → datasets
        n_obs = Column(Integer)  # total cells/spots at this phase
        n_vars = Column(Integer)  # genes/proteins at this phase
        n_samples = Column(Integer)  # number of samples at this phase
        notes = Column(Text)  # free-text notes
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        project = relationship("Project", back_populates="phases")

    class DataSource(Base):
        """A raw data source — HPC directory, public dataset, or external repository.

        Distinct from :class:`Project` (lab work) and :class:`Dataset` (processed
        checkpoints). Projects reference DataSources via the
        :class:`ProjectDataSource` join table.
        """

        __tablename__ = "data_sources"

        id = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String, unique=True, nullable=False)
        description = Column(Text)
        uri = Column(Text, nullable=False)  # HPC path, URL, GEO accession, DOI, etc.
        platform = Column(String)  # cosmx, xenium, imc, visium_hd, scrna, ...
        domain = Column(String)  # spatial_transcriptomics | spatial_proteomics | ...
        imaging_modality = Column(String)  # sequencing_based | probe_based | ...
        source_type = Column(
            String
        )  # hpc_lab | hpc_collaborator | public_10x | public_geo | public_zenodo | public_portal
        organism = Column(String)  # human | mouse | ...
        tissue = Column(String)  # colon | brain | ...
        disease = Column(String)  # IBD | breast_cancer | ...
        n_samples = Column(Integer)
        n_cells = Column(Integer)
        publication = Column(String)  # DOI or citation string
        access_notes = Column(Text)  # e.g. "cayuga:/athena/project-saha/; contact jip2007"
        status = Column(
            String, default="available"
        )  # available | restricted | pending_download | archived
        created_at = Column(String, default=_utcnow)

        project_links = relationship(
            "ProjectDataSource", back_populates="data_source", cascade="all, delete-orphan"
        )

    class ProjectDataSource(Base):
        """Many-to-many join between Projects and DataSources."""

        __tablename__ = "project_data_sources"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
        data_source_id = Column(
            Integer, ForeignKey("data_sources.id", ondelete="CASCADE"), nullable=False
        )
        role = Column(String, default="input")  # input | reference | supplementary
        notes = Column(Text)

        project = relationship("Project", back_populates="data_source_links")
        data_source = relationship("DataSource", back_populates="project_links")

    # Wire back-ref on Project
    Project.data_source_links = relationship(
        "ProjectDataSource", back_populates="project", cascade="all, delete-orphan"
    )

    # ------------------------------------------------------------------
    # BioData hierarchy (migration 0005)
    # ------------------------------------------------------------------

    class Subject(Base):
        """Cross-project de-identified subject/patient record."""

        __tablename__ = "subjects"

        id = Column(Integer, primary_key=True, autoincrement=True)
        subject_id = Column(String, unique=True, nullable=False)
        organism = Column(String, default="human")
        # Standardized clinical columns
        sex = Column(String)  # M, F, unknown
        age_at_collection = Column(Float)
        diagnosis = Column(String)
        diagnosis_code = Column(String)  # ICD-10 or OncoTree
        disease_stage = Column(String)
        treatment_status = Column(String)  # treatment_naive | on_treatment | post_treatment
        tissue_of_origin = Column(String)
        cause_of_death = Column(String)
        survival_days = Column(Float)
        vital_status = Column(String)  # alive | deceased | unknown
        # Flexible overflow
        clinical_metadata_json = Column(Text)
        created_at = Column(String, default=_utcnow)

        samples = relationship("Sample", back_populates="subject", cascade="all, delete-orphan")
        project_links = relationship(
            "SubjectProjectLink", back_populates="subject", cascade="all, delete-orphan"
        )

    class Sample(Base):
        """Physical specimen linking a subject to data within a project."""

        __tablename__ = "samples"

        id = Column(Integer, primary_key=True, autoincrement=True)
        sample_id = Column(String, nullable=False)
        subject_id = Column(Integer, ForeignKey("subjects.id", ondelete="SET NULL"), nullable=True)
        project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=True)
        tissue = Column(String)
        tissue_region = Column(String)
        collection_date = Column(String)
        fixation_method = Column(String)  # FFPE | fresh_frozen | OCT | PFA
        sample_type = Column(String)  # biopsy | resection | TMA | organoid | xenograft
        batch = Column(String)
        notes = Column(Text)
        created_at = Column(String, default=_utcnow)

        subject = relationship("Subject", back_populates="samples")
        project = relationship("Project")

    class SubjectProjectLink(Base):
        """Many-to-many join between Subjects and Projects."""

        __tablename__ = "subject_project_links"

        id = Column(Integer, primary_key=True, autoincrement=True)
        subject_id = Column(Integer, ForeignKey("subjects.id", ondelete="CASCADE"), nullable=False)
        project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
        role = Column(String, default="enrolled")  # enrolled | reference | control
        notes = Column(Text)

        __table_args__ = (UniqueConstraint("subject_id", "project_id", name="uq_subject_project"),)

        subject = relationship("Subject", back_populates="project_links")
        project = relationship("Project")

    class BioData(Base):
        """Joined Table Inheritance base for typed biological data objects."""

        __tablename__ = "bio_data"

        id = Column(Integer, primary_key=True, autoincrement=True)
        type = Column(String, nullable=False)  # JTI discriminator
        project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
        subject_id = Column(Integer, ForeignKey("subjects.id", ondelete="SET NULL"), nullable=True)
        sample_id = Column(Integer, ForeignKey("samples.id", ondelete="SET NULL"), nullable=True)
        # Classification
        category = Column(String)
        subcategory = Column(String)
        platform = Column(String)
        modality = Column(String)
        measurement = Column(String)
        resolution = Column(String)
        spatial = Column(Boolean)
        # Data tracking (replaces Dataset for new data)
        uri = Column(Text, nullable=False)
        format = Column(String)
        status = Column(String, default="pending")
        file_role = Column(String, default="primary")
        validated = Column(Boolean, default=False)
        n_obs = Column(Integer)
        n_vars = Column(Integer)
        phase = Column(String)
        size_mb = Column(Float)
        md5 = Column(String)
        # Provenance
        legacy_dataset_id = Column(
            Integer, ForeignKey("datasets.id", ondelete="SET NULL"), nullable=True
        )
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        project = relationship("Project")

        __mapper_args__ = {
            "polymorphic_on": type,
            "polymorphic_identity": "base",
        }

    class BioImage(BioData):
        """Image-specific BioData (IMC, CODEX, H&E, IF, etc.)."""

        __tablename__ = "bio_images"

        id = Column(Integer, ForeignKey("bio_data.id", ondelete="CASCADE"), primary_key=True)
        image_type = Column(String)  # he, fluorescence, multiplexed, brightfield, phase_contrast
        n_channels = Column(Integer)
        pixel_size_um = Column(Float)
        width_px = Column(Integer)
        height_px = Column(Integer)
        staining_protocol = Column(String)  # H&E, IHC, IF, CODEX, IMC
        channel_names_json = Column(Text)  # JSON list of channel/marker names

        __mapper_args__ = {"polymorphic_identity": "image"}

    class RNASeqData(BioData):
        """RNA-seq-specific BioData (bulk and single-cell)."""

        __tablename__ = "rnaseq_data"

        id = Column(Integer, ForeignKey("bio_data.id", ondelete="CASCADE"), primary_key=True)
        library_type = Column(String)  # bulk | single_cell
        sequencing_platform = Column(String)  # illumina_novaseq, bgi_dnbseq, etc.
        chemistry = Column(String)  # chromium_v3, evercode_v3, etc.
        read_length = Column(Integer)
        n_reads = Column(Integer)
        reference_genome = Column(String)  # GRCh38, GRCm39, etc.
        gene_annotation = Column(String)  # gencode_v44, ensembl_110, etc.

        __mapper_args__ = {"polymorphic_identity": "rnaseq"}

    class SpatialSeqData(BioData):
        """Spatial sequencing-specific BioData (Visium, Xenium, CosMx, etc.)."""

        __tablename__ = "spatial_seq_data"

        id = Column(Integer, ForeignKey("bio_data.id", ondelete="CASCADE"), primary_key=True)
        spatial_resolution = Column(String)  # spot | single_cell
        panel_size = Column(Integer)  # number of targeted genes
        bin_size_um = Column(Float)  # for Visium HD bin-level
        fov_count = Column(Integer)  # for CosMx
        coordinate_system = Column(String)  # pixel | micron
        tissue_area_mm2 = Column(Float)
        sequencing_platform = Column(String)

        __mapper_args__ = {"polymorphic_identity": "spatial_seq"}

    class EpigenomicsData(BioData):
        """Epigenomics-specific BioData (ATAC-seq, ChIP-seq, CUT&Tag, etc.)."""

        __tablename__ = "epigenomics_data"

        id = Column(Integer, ForeignKey("bio_data.id", ondelete="CASCADE"), primary_key=True)
        assay_type = Column(String)  # atac_seq, chip_seq, cut_and_tag, methylation, etc.
        target_protein = Column(String)  # for ChIP/CUT&Tag
        n_peaks = Column(Integer)
        n_fragments = Column(Integer)
        genome_coverage = Column(Float)

        __mapper_args__ = {"polymorphic_identity": "epigenomics"}

    class GenomeSeqData(BioData):
        """Genome sequencing-specific BioData (WGS, WES, targeted panels)."""

        __tablename__ = "genome_seq_data"

        id = Column(Integer, ForeignKey("bio_data.id", ondelete="CASCADE"), primary_key=True)
        sequencing_type = Column(String)  # wgs | wes | targeted_panel | amplicon
        panel_name = Column(String)
        coverage_mean = Column(Float)
        reference_genome = Column(String)  # GRCh38, etc.
        sequencing_platform = Column(String)  # illumina_novaseq, ont_promethion, etc.

        __mapper_args__ = {"polymorphic_identity": "genome_seq"}

    return (
        Project,
        Dataset,
        SlurmJob,
        AgentTask,
        ProjectPhase,
        DataSource,
        ProjectDataSource,
        Subject,
        Sample,
        SubjectProjectLink,
        BioData,
        BioImage,
        RNASeqData,
        SpatialSeqData,
        EpigenomicsData,
        GenomeSeqData,
    )


# ---------------------------------------------------------------------------
# Registry class
# ---------------------------------------------------------------------------


class Registry:
    """Project and data registry backed by SQLite or PostgreSQL.

    Parameters
    ----------
    db_url
        SQLAlchemy DB URL.  Defaults to env var ``SC_TOOLS_REGISTRY_URL``
        or ``sqlite:///~/.sc_tools/registry.db``.
    """

    def __init__(self, db_url: str | None = None) -> None:
        sa = _require_sqlalchemy()
        from sqlalchemy.orm import DeclarativeBase

        class Base(DeclarativeBase):
            pass

        self._Base = Base
        (
            self._Project,
            self._Dataset,
            self._SlurmJob,
            self._AgentTask,
            self._ProjectPhase,
            self._DataSource,
            self._ProjectDataSource,
            self._Subject,
            self._Sample,
            self._SubjectProjectLink,
            self._BioData,
            self._BioImage,
            self._RNASeqData,
            self._SpatialSeqData,
            self._EpigenomicsData,
            self._GenomeSeqData,
        ) = _build_models(Base)

        url = db_url or _get_db_url()
        self._engine = sa.create_engine(url, echo=False)
        Base.metadata.create_all(self._engine)
        logger.info("Registry connected: %s", url)

    def _session(self) -> Any:
        from sqlalchemy.orm import Session

        return Session(self._engine)

    # ------------------------------------------------------------------
    # Projects
    # ------------------------------------------------------------------

    def add_project(
        self,
        name: str,
        platform: str | None = None,
        data_type: str | None = None,
        *,
        domain: str | None = None,
        imaging_modality: str | None = None,
        project_type: str = "internal",
        visibility: str = "private",
    ) -> int:
        """Register a new project. Returns project id.

        If a project with the same name exists its id is returned without
        creating a duplicate.

        Parameters
        ----------
        name:
            Unique project name (e.g. ``"ggo_visium"``).
        platform:
            Technology platform string (e.g. ``"visium"``, ``"imc"``).
        data_type:
            Data type (may match platform for single-modality projects).
        domain:
            High-level domain: ``spatial_transcriptomics`` |
            ``spatial_proteomics`` | ``imaging`` | ``single_cell`` | ``bulk``.
        imaging_modality:
            Imaging modality: ``brightfield`` | ``fluorescence`` |
            ``multiplexed_fluorescence`` | ``probe_based`` |
            ``mass_spec_imaging`` | ``sequencing_based``.
        project_type:
            ``"internal"`` (lab-led) or ``"external"`` (collaboration / contract).
        visibility:
            ``"private"`` (restricted access) or ``"public"`` (shareable / published).
        """
        with self._session() as sess:
            existing = sess.query(self._Project).filter_by(name=name).first()
            if existing:
                logger.info("Project '%s' already registered (id=%d)", name, existing.id)
                return existing.id
            proj = self._Project(
                name=name,
                platform=platform,
                data_type=data_type,
                domain=domain,
                imaging_modality=imaging_modality,
                project_type=project_type,
                visibility=visibility,
            )
            sess.add(proj)
            sess.commit()
            sess.refresh(proj)
            logger.info("Registered project '%s' (id=%d)", name, proj.id)
            return proj.id

    def get_project(self, name: str) -> dict[str, Any] | None:
        """Return project dict or None if not found."""
        with self._session() as sess:
            row = sess.query(self._Project).filter_by(name=name).first()
            return self._to_dict(row) if row else None

    def list_projects(self) -> list[dict[str, Any]]:
        """Return all registered projects."""
        with self._session() as sess:
            return [self._to_dict(r) for r in sess.query(self._Project).all()]

    def mark_phase_complete(self, project_name: str, phase: str) -> None:
        """Append *phase* to the project's completed-phases list and upsert project_phases."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")
            phases: list[str] = json.loads(proj.phases_complete or "[]")
            if phase not in phases:
                phases.append(phase)
                proj.phases_complete = json.dumps(phases)
                sess.commit()
        # Also update the project_phases table
        self.upsert_phase(project_name, phase, status="complete")

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------

    def register_dataset(
        self,
        project_name: str,
        phase: str,
        uri: str,
        *,
        sample_id: str | None = None,
        fmt: str = "h5ad",
        size_mb: float | None = None,
        md5: str | None = None,
        status: str = "ready",
        file_role: str = "primary",
        validated: bool = False,
        n_obs: int | None = None,
        n_vars: int | None = None,
    ) -> int:
        """Register an AnnData checkpoint. Returns dataset id.

        .. deprecated::
            Use :meth:`register_biodata` for new data. This method creates
            a legacy ``datasets`` row. In a future release it will internally
            create a BioData row instead.

        Parameters
        ----------
        project_name:
            Project name (must already exist).
        phase:
            Phase slug (e.g. ``"qc_filter"``, ``"scoring"``).
        uri:
            Path or URI to the checkpoint file.
        sample_id:
            Optional sample identifier (for per-sample ingest_load checkpoints).
        fmt:
            File format: ``"h5ad"`` (default), ``"zarr"``, ``"tiff"``, ``"tsv"``.
        file_role:
            Role of this file: ``"primary"`` | ``"supplementary"`` |
            ``"entry_point"`` | ``"spatialdata"`` | ``"image"`` | ``"metadata"``.
        validated:
            True if ``validate_checkpoint()`` has passed for this file.
        n_obs:
            Number of cells/spots in this dataset (for reporting).
        n_vars:
            Number of genes/proteins in this dataset (for reporting).
        """
        import warnings

        warnings.warn(
            "register_dataset() is deprecated; use register_biodata() for new data. "
            "register_dataset() will internally create BioData rows in a future release.",
            DeprecationWarning,
            stacklevel=2,
        )
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found. Call add_project() first.")
            ds = self._Dataset(
                project_id=proj.id,
                sample_id=sample_id,
                phase=phase,
                uri=uri,
                format=fmt,
                size_mb=size_mb,
                md5=md5,
                status=status,
                file_role=file_role,
                validated=validated,
                n_obs=n_obs,
                n_vars=n_vars,
            )
            sess.add(ds)
            sess.commit()
            sess.refresh(ds)
            ds_id = ds.id
            proj_platform = proj.platform
            logger.info("Registered dataset %s/%s at %s (id=%d)", project_name, phase, uri, ds_id)

        # Phase B dual-write: also create a BioData row
        try:
            from sc_tools.biodata import platform_for_project

            pspec = platform_for_project(proj_platform) if proj_platform else None
            inferred_category = pspec.biodata_type if pspec else "spatial_seq"
            bio_id = self.register_biodata(
                project_name=project_name,
                category=inferred_category,
                platform=proj_platform or "unknown",
                uri=uri,
                fmt=fmt,
                phase=phase,
                file_role=file_role,
                n_obs=n_obs,
                n_vars=n_vars,
                status=status,
            )
            # Link forward: dataset -> biodata
            with self._session() as sess:
                ds_row = sess.get(self._Dataset, ds_id)
                if ds_row is not None:
                    ds_row.bio_data_id = bio_id
                    sess.commit()
        except Exception:
            logger.warning("Dual-write to BioData failed for dataset %d", ds_id)

        return ds_id

    def get_dataset_uri(
        self,
        project_name: str,
        phase: str,
        sample_id: str | None = None,
    ) -> str | None:
        """Return the URI for a checkpoint (most recently registered).

        Returns None if no matching dataset is found.
        """
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return None
            q = sess.query(self._Dataset).filter_by(project_id=proj.id, phase=phase)
            if sample_id is not None:
                q = q.filter_by(sample_id=sample_id)
            row = q.order_by(self._Dataset.id.desc()).first()
            return row.uri if row else None

    def list_datasets(
        self,
        project_name: str | None = None,
        phase: str | None = None,
    ) -> list[dict[str, Any]]:
        """List datasets, optionally filtered by project and/or phase."""
        with self._session() as sess:
            q = sess.query(self._Dataset)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj:
                    q = q.filter_by(project_id=proj.id)
                else:
                    return []
            if phase is not None:
                q = q.filter_by(phase=phase)
            return [self._to_dict(r) for r in q.all()]

    def update_dataset_status(self, dataset_id: int, status: str) -> None:
        """Update dataset status (pending, ready, archived, error)."""
        with self._session() as sess:
            ds = sess.get(self._Dataset, dataset_id)
            if ds is None:
                raise ValueError(f"Dataset id={dataset_id} not found")
            ds.status = status
            ds.updated_at = _utcnow()
            sess.commit()

    def mark_dataset_validated(self, dataset_id: int) -> None:
        """Set ``datasets.validated = True`` for the given dataset id."""
        with self._session() as sess:
            ds = sess.get(self._Dataset, dataset_id)
            if ds is None:
                raise ValueError(f"Dataset id={dataset_id} not found")
            ds.validated = True
            ds.updated_at = _utcnow()
            sess.commit()

    # ------------------------------------------------------------------
    # Project phases
    # ------------------------------------------------------------------

    def upsert_phase(
        self,
        project_name: str,
        phase: str,
        *,
        status: str = "in_progress",
        entry_phase: bool = False,
        primary_dataset_id: int | None = None,
        n_obs: int | None = None,
        n_vars: int | None = None,
        n_samples: int | None = None,
        notes: str | None = None,
    ) -> None:
        """Insert or update a ``project_phases`` row (upsert on composite PK).

        Parameters
        ----------
        project_name:
            Project name (must already exist).
        phase:
            Phase slug (e.g. ``"qc_filter"``).
        status:
            Pipeline status: ``"not_started"`` | ``"in_progress"`` |
            ``"complete"`` | ``"failed"`` | ``"skipped"``.
        entry_phase:
            True if this is the phase at which the pipeline was loaded
            (e.g. collaborator provided a pre-celltyped AnnData).
        primary_dataset_id:
            FK to ``datasets.id`` for the primary checkpoint of this phase.
        n_obs:
            Total cells/spots at this phase.
        n_vars:
            Genes/proteins at this phase.
        n_samples:
            Number of samples at this phase.
        notes:
            Free-text notes (e.g. ``"9-method benchmark; Z-score+Harmony selected"``).
        """
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")

            existing = (
                sess.query(self._ProjectPhase).filter_by(project_id=proj.id, phase=phase).first()
            )
            now = _utcnow()
            if existing is None:
                row = self._ProjectPhase(
                    project_id=proj.id,
                    phase=phase,
                    status=status,
                    entry_phase=entry_phase,
                    primary_dataset_id=primary_dataset_id,
                    n_obs=n_obs,
                    n_vars=n_vars,
                    n_samples=n_samples,
                    notes=notes,
                    created_at=now,
                    updated_at=now,
                )
                sess.add(row)
            else:
                existing.status = status
                existing.updated_at = now
                if entry_phase:
                    existing.entry_phase = entry_phase
                if primary_dataset_id is not None:
                    existing.primary_dataset_id = primary_dataset_id
                if n_obs is not None:
                    existing.n_obs = n_obs
                if n_vars is not None:
                    existing.n_vars = n_vars
                if n_samples is not None:
                    existing.n_samples = n_samples
                if notes is not None:
                    existing.notes = notes
            sess.commit()

    def get_phase(self, project_name: str, phase: str) -> dict[str, Any] | None:
        """Return the ``project_phases`` row for ``(project, phase)`` as dict, or None."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return None
            row = sess.query(self._ProjectPhase).filter_by(project_id=proj.id, phase=phase).first()
            return self._to_dict(row) if row else None

    def list_phases(self, project_name: str) -> list[dict[str, Any]]:
        """Return all ``project_phases`` rows for a project, ordered by phase slug."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return []
            rows = (
                sess.query(self._ProjectPhase)
                .filter_by(project_id=proj.id)
                .order_by(self._ProjectPhase.phase)
                .all()
            )
            return [self._to_dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Data sources
    # ------------------------------------------------------------------

    def register_data_source(
        self,
        name: str,
        uri: str,
        *,
        description: str | None = None,
        platform: str | None = None,
        domain: str | None = None,
        imaging_modality: str | None = None,
        source_type: str | None = None,
        organism: str | None = None,
        tissue: str | None = None,
        disease: str | None = None,
        n_samples: int | None = None,
        n_cells: int | None = None,
        publication: str | None = None,
        access_notes: str | None = None,
        status: str = "available",
    ) -> int:
        """Register a data source (HPC directory, public dataset, external repo).

        Returns the data source id. Idempotent: if a source with the same
        ``name`` already exists, its id is returned without creating a duplicate.

        Parameters
        ----------
        name:
            Unique identifier (e.g. ``"saha_ibd_cosmx"``).
        uri:
            Primary location — HPC path (``"cayuga:/athena/project-saha/data_IBD"``),
            URL, GEO accession, DOI, etc.
        source_type:
            ``"hpc_lab"`` | ``"hpc_collaborator"`` | ``"public_10x"`` |
            ``"public_geo"`` | ``"public_zenodo"`` | ``"public_portal"``.
        status:
            ``"available"`` | ``"restricted"`` | ``"pending_download"`` | ``"archived"``.
        """
        with self._session() as sess:
            existing = sess.query(self._DataSource).filter_by(name=name).first()
            if existing:
                logger.info("DataSource '%s' already registered (id=%d)", name, existing.id)
                return existing.id
            src = self._DataSource(
                name=name,
                uri=uri,
                description=description,
                platform=platform,
                domain=domain,
                imaging_modality=imaging_modality,
                source_type=source_type,
                organism=organism,
                tissue=tissue,
                disease=disease,
                n_samples=n_samples,
                n_cells=n_cells,
                publication=publication,
                access_notes=access_notes,
                status=status,
            )
            sess.add(src)
            sess.commit()
            sess.refresh(src)
            logger.info("Registered data source '%s' (id=%d)", name, src.id)
            return src.id

    def get_data_source(self, name: str) -> dict[str, Any] | None:
        """Return data source dict or None if not found."""
        with self._session() as sess:
            row = sess.query(self._DataSource).filter_by(name=name).first()
            return self._to_dict(row) if row else None

    def list_data_sources(
        self,
        platform: str | None = None,
        source_type: str | None = None,
        disease: str | None = None,
        tissue: str | None = None,
    ) -> list[dict[str, Any]]:
        """Return data sources, optionally filtered."""
        with self._session() as sess:
            q = sess.query(self._DataSource)
            if platform is not None:
                q = q.filter_by(platform=platform)
            if source_type is not None:
                q = q.filter_by(source_type=source_type)
            if disease is not None:
                q = q.filter(self._DataSource.disease.ilike(f"%{disease}%"))
            if tissue is not None:
                q = q.filter(self._DataSource.tissue.ilike(f"%{tissue}%"))
            return [self._to_dict(r) for r in q.order_by(self._DataSource.name).all()]

    def link_project_data_source(
        self,
        project_name: str,
        data_source_name: str,
        *,
        role: str = "input",
        notes: str | None = None,
    ) -> int:
        """Link a project to a data source. Returns link id.

        Parameters
        ----------
        role:
            ``"input"`` (primary data used in this project) |
            ``"reference"`` (used as reference, e.g. scRNA ref for deconvolution) |
            ``"supplementary"`` (additional context).
        """
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found.")
            src = sess.query(self._DataSource).filter_by(name=data_source_name).first()
            if src is None:
                raise ValueError(f"DataSource '{data_source_name}' not found.")
            # Idempotent — skip if link already exists
            existing = (
                sess.query(self._ProjectDataSource)
                .filter_by(project_id=proj.id, data_source_id=src.id)
                .first()
            )
            if existing:
                return existing.id
            link = self._ProjectDataSource(
                project_id=proj.id,
                data_source_id=src.id,
                role=role,
                notes=notes,
            )
            sess.add(link)
            sess.commit()
            sess.refresh(link)
            logger.info(
                "Linked project '%s' → data source '%s' (role=%s)",
                project_name,
                data_source_name,
                role,
            )
            return link.id

    def list_project_data_sources(self, project_name: str) -> list[dict[str, Any]]:
        """Return all data sources linked to a project."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return []
            links = sess.query(self._ProjectDataSource).filter_by(project_id=proj.id).all()
            result = []
            for lnk in links:
                src = sess.get(self._DataSource, lnk.data_source_id)
                if src:
                    d = self._to_dict(src)
                    d["link_role"] = lnk.role
                    d["link_notes"] = lnk.notes
                    result.append(d)
            return result

    def delete_project(self, name: str) -> bool:
        """Delete a project and all its cascade rows. Returns True if deleted."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=name).first()
            if proj is None:
                return False
            sess.delete(proj)
            sess.commit()
            logger.info("Deleted project '%s'", name)
            return True

    # ------------------------------------------------------------------
    # SLURM jobs
    # ------------------------------------------------------------------

    def register_job(
        self,
        project_name: str,
        sample_id: str | None,
        phase: str,
        cluster: str,
        slurm_job_id: str,
    ) -> int:
        """Record a submitted SLURM job. Returns job id."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found.")
            job = self._SlurmJob(
                project_id=proj.id,
                sample_id=sample_id,
                phase=phase,
                cluster=cluster,
                slurm_job_id=slurm_job_id,
                status="submitted",
                submitted_at=_utcnow(),
            )
            sess.add(job)
            sess.commit()
            sess.refresh(job)
            logger.info("Registered SLURM job %s on %s (id=%d)", slurm_job_id, cluster, job.id)
            return job.id

    def update_job_status(
        self,
        slurm_job_id: str,
        status: str,
        *,
        error: str | None = None,
        log_uri: str | None = None,
    ) -> None:
        """Update status for a SLURM job identified by its slurm_job_id."""
        with self._session() as sess:
            job = sess.query(self._SlurmJob).filter_by(slurm_job_id=slurm_job_id).first()
            if job is None:
                raise ValueError(f"SLURM job '{slurm_job_id}' not found in registry")
            job.status = status
            if status in ("completed", "failed"):
                job.finished_at = _utcnow()
            if error is not None:
                job.error_msg = error
            if log_uri is not None:
                job.log_uri = log_uri
            sess.commit()

    def list_active_jobs(self) -> list[dict[str, Any]]:
        """Return all jobs with status submitted or running."""
        with self._session() as sess:
            rows = (
                sess.query(self._SlurmJob)
                .filter(self._SlurmJob.status.in_(["submitted", "running"]))
                .all()
            )
            return [self._to_dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Agent tasks
    # ------------------------------------------------------------------

    def start_task(
        self,
        task_type: str,
        project_name: str,
        inputs: dict[str, Any] | None = None,
    ) -> int:
        """Record the start of an agent task. Returns task id."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found.")
            task = self._AgentTask(
                task_type=task_type,
                project_id=proj.id,
                status="running",
                inputs_json=json.dumps(inputs or {}),
                started_at=_utcnow(),
            )
            sess.add(task)
            sess.commit()
            sess.refresh(task)
            return task.id

    def finish_task(
        self,
        task_id: int,
        *,
        outputs: dict[str, Any] | None = None,
        error: str | None = None,
    ) -> None:
        """Mark an agent task as completed or failed."""
        with self._session() as sess:
            task = sess.get(self._AgentTask, task_id)
            if task is None:
                raise ValueError(f"AgentTask id={task_id} not found")
            task.status = "failed" if error else "completed"
            task.finished_at = _utcnow()
            if outputs is not None:
                task.outputs_json = json.dumps(outputs)
            if error is not None:
                task.error = error
            sess.commit()

    def list_running_tasks(self) -> list[dict[str, Any]]:
        """Return all running agent tasks."""
        with self._session() as sess:
            rows = sess.query(self._AgentTask).filter_by(status="running").all()
            return [self._to_dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Subjects
    # ------------------------------------------------------------------

    def add_subject(
        self,
        subject_id: str,
        organism: str = "human",
        **clinical_kwargs: Any,
    ) -> int:
        """Register a de-identified subject. Returns subject DB id.

        Idempotent: if a subject with the same ``subject_id`` exists, its
        DB id is returned without creating a duplicate.

        Parameters
        ----------
        subject_id:
            Unique de-identified identifier (e.g. ``"PT001"``).
        organism:
            ``"human"`` (default), ``"mouse"``, ``"rat"``, etc.
        **clinical_kwargs:
            Any standardized clinical column (``sex``, ``age_at_collection``,
            ``diagnosis``, ``diagnosis_code``, ``disease_stage``,
            ``treatment_status``, ``tissue_of_origin``, ``cause_of_death``,
            ``survival_days``, ``vital_status``) or ``clinical_metadata_json``
            (JSON string for overflow fields).
        """
        _standard_cols = {
            "sex",
            "age_at_collection",
            "diagnosis",
            "diagnosis_code",
            "disease_stage",
            "treatment_status",
            "tissue_of_origin",
            "cause_of_death",
            "survival_days",
            "vital_status",
            "clinical_metadata_json",
        }
        with self._session() as sess:
            existing = sess.query(self._Subject).filter_by(subject_id=subject_id).first()
            if existing:
                logger.info("Subject '%s' already registered (id=%d)", subject_id, existing.id)
                return existing.id
            kwargs_filtered = {k: v for k, v in clinical_kwargs.items() if k in _standard_cols}
            subj = self._Subject(
                subject_id=subject_id,
                organism=organism,
                **kwargs_filtered,
            )
            sess.add(subj)
            sess.commit()
            sess.refresh(subj)
            logger.info("Registered subject '%s' (id=%d)", subject_id, subj.id)
            return subj.id

    def get_subject(self, subject_id: str) -> dict[str, Any] | None:
        """Return subject dict by de-identified ID, or None."""
        with self._session() as sess:
            row = sess.query(self._Subject).filter_by(subject_id=subject_id).first()
            return self._to_dict(row) if row else None

    def list_subjects(
        self,
        project_name: str | None = None,
        diagnosis: str | None = None,
        tissue: str | None = None,
    ) -> list[dict[str, Any]]:
        """Return subjects, optionally filtered."""
        with self._session() as sess:
            q = sess.query(self._Subject)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                link_subquery = (
                    sess.query(self._SubjectProjectLink.subject_id)
                    .filter_by(project_id=proj.id)
                    .scalar_subquery()
                )
                q = q.filter(self._Subject.id.in_(link_subquery))
            if diagnosis is not None:
                q = q.filter(self._Subject.diagnosis.ilike(f"%{diagnosis}%"))
            if tissue is not None:
                q = q.filter(self._Subject.tissue_of_origin.ilike(f"%{tissue}%"))
            return [self._to_dict(r) for r in q.order_by(self._Subject.subject_id).all()]

    def link_subject_to_project(
        self,
        subject_id: str,
        project_name: str,
        role: str = "enrolled",
    ) -> int:
        """Link a subject to a project. Returns link DB id.

        Idempotent: returns existing link id if already linked.
        """
        with self._session() as sess:
            subj = sess.query(self._Subject).filter_by(subject_id=subject_id).first()
            if subj is None:
                raise ValueError(f"Subject '{subject_id}' not found in registry")
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")
            existing = (
                sess.query(self._SubjectProjectLink)
                .filter_by(subject_id=subj.id, project_id=proj.id)
                .first()
            )
            if existing:
                return existing.id
            link = self._SubjectProjectLink(
                subject_id=subj.id,
                project_id=proj.id,
                role=role,
            )
            sess.add(link)
            sess.commit()
            sess.refresh(link)
            return link.id

    # ------------------------------------------------------------------
    # Samples
    # ------------------------------------------------------------------

    def add_sample(
        self,
        sample_id: str,
        subject_id: str | None = None,
        project_name: str | None = None,
        *,
        tissue: str | None = None,
        tissue_region: str | None = None,
        collection_date: str | None = None,
        fixation_method: str | None = None,
        sample_type: str | None = None,
        batch: str | None = None,
        notes: str | None = None,
    ) -> int:
        """Register a sample. Returns sample DB id.

        Parameters
        ----------
        sample_id:
            Sample identifier (e.g. ``"S001_A1"``).
        subject_id:
            De-identified subject ID (optional). Must already exist.
        project_name:
            Project name (optional). Must already exist.
        """
        with self._session() as sess:
            subj_db_id = None
            if subject_id is not None:
                subj = sess.query(self._Subject).filter_by(subject_id=subject_id).first()
                if subj is None:
                    raise ValueError(f"Subject '{subject_id}' not found in registry")
                subj_db_id = subj.id
            proj_db_id = None
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    raise ValueError(f"Project '{project_name}' not found in registry")
                proj_db_id = proj.id
            sample = self._Sample(
                sample_id=sample_id,
                subject_id=subj_db_id,
                project_id=proj_db_id,
                tissue=tissue,
                tissue_region=tissue_region,
                collection_date=collection_date,
                fixation_method=fixation_method,
                sample_type=sample_type,
                batch=batch,
                notes=notes,
            )
            sess.add(sample)
            sess.commit()
            sess.refresh(sample)
            logger.info("Registered sample '%s' (id=%d)", sample_id, sample.id)
            return sample.id

    def get_sample(
        self,
        sample_id: str,
        project_name: str | None = None,
    ) -> dict[str, Any] | None:
        """Return sample dict by sample_id, optionally scoped to project."""
        with self._session() as sess:
            q = sess.query(self._Sample).filter_by(sample_id=sample_id)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return None
                q = q.filter_by(project_id=proj.id)
            row = q.first()
            return self._to_dict(row) if row else None

    def list_samples(
        self,
        project_name: str | None = None,
        subject_id: str | None = None,
        batch: str | None = None,
    ) -> list[dict[str, Any]]:
        """Return samples, optionally filtered."""
        with self._session() as sess:
            q = sess.query(self._Sample)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                q = q.filter_by(project_id=proj.id)
            if subject_id is not None:
                subj = sess.query(self._Subject).filter_by(subject_id=subject_id).first()
                if subj is None:
                    return []
                q = q.filter_by(subject_id=subj.id)
            if batch is not None:
                q = q.filter_by(batch=batch)
            return [self._to_dict(r) for r in q.order_by(self._Sample.sample_id).all()]

    # ------------------------------------------------------------------
    # BioData (polymorphic)
    # ------------------------------------------------------------------

    _BIODATA_TYPE_MAP: dict[str, str] = {
        "image": "_BioImage",
        "rnaseq": "_RNASeqData",
        "spatial_seq": "_SpatialSeqData",
        "epigenomics": "_EpigenomicsData",
        "genome_seq": "_GenomeSeqData",
    }

    def register_biodata(
        self,
        project_name: str,
        category: str,
        platform: str,
        uri: str,
        *,
        subcategory: str | None = None,
        measurement: str | None = None,
        resolution: str | None = None,
        spatial: bool | None = None,
        fmt: str | None = None,
        status: str = "pending",
        file_role: str = "primary",
        validated: bool = False,
        n_obs: int | None = None,
        n_vars: int | None = None,
        phase: str | None = None,
        size_mb: float | None = None,
        md5: str | None = None,
        subject_id: str | None = None,
        sample_db_id: int | None = None,
        **type_kwargs: Any,
    ) -> int:
        """Register a typed BioData object. Returns bio_data id.

        The ``category`` determines the JTI subtype: ``"image"`` creates a
        :class:`BioImage`, ``"spatial_seq"`` creates :class:`SpatialSeqData`,
        etc. Extra ``**type_kwargs`` are passed to the subtype constructor
        (e.g. ``image_type="multiplexed"`` for BioImage).

        If ``platform`` is in ``KNOWN_PLATFORMS``, defaults for ``subcategory``,
        ``measurement``, ``resolution``, ``spatial`` are auto-filled (but
        explicit values take precedence).
        """
        # Auto-fill from platform registry
        modality: str | None = None
        pspec = None
        try:
            from sc_tools.biodata import get_platform

            pspec = get_platform(platform)
            if subcategory is None:
                subcategory = pspec.subcategory
            if measurement is None:
                measurement = pspec.measurement
            if resolution is None:
                resolution = pspec.resolution
            if spatial is None:
                spatial = pspec.spatial
            if category is None:
                category = pspec.category
            modality = pspec.modality or None
            # Auto-fill child columns from PlatformSpec.defaults
            if pspec.defaults:
                for k, v in pspec.defaults.items():
                    if k not in type_kwargs:
                        type_kwargs[k] = v
        except (KeyError, ImportError):
            pass

        # Auto-fill format from URI extension
        if fmt is None and uri:
            ext = uri.rsplit(".", 1)[-1].lower() if "." in uri else ""
            if ext in ("h5ad", "zarr", "tiff", "tsv", "fastq", "bam", "bed", "csv"):
                fmt = ext

        # Determine model class
        model_attr = self._BIODATA_TYPE_MAP.get(category)
        if model_attr is None:
            raise ValueError(
                f"Unknown BioData category '{category}'. "
                f"Must be one of: {list(self._BIODATA_TYPE_MAP.keys())}"
            )
        ModelClass = getattr(self, model_attr)

        # Resolve FKs
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found. Call add_project() first.")
            subj_db_id = None
            if subject_id is not None:
                subj = sess.query(self._Subject).filter_by(subject_id=subject_id).first()
                if subj is not None:
                    subj_db_id = subj.id

            # Filter type_kwargs to only valid columns for the model
            valid_cols = {c.name for c in ModelClass.__table__.columns} - {
                "id",
                "type",
                "project_id",
                "subject_id",
                "sample_id",
                "category",
                "subcategory",
                "platform",
                "modality",
                "measurement",
                "resolution",
                "spatial",
                "uri",
                "format",
                "status",
                "file_role",
                "validated",
                "n_obs",
                "n_vars",
                "phase",
                "size_mb",
                "md5",
                "legacy_dataset_id",
                "created_at",
                "updated_at",
            }
            filtered_kwargs = {k: v for k, v in type_kwargs.items() if k in valid_cols}
            unknown_kwargs = set(type_kwargs.keys()) - valid_cols
            if unknown_kwargs:
                logger.warning(
                    "register_biodata: ignoring unknown kwargs for %s: %s. "
                    "Valid type-specific columns: %s",
                    category,
                    unknown_kwargs,
                    valid_cols,
                )

            obj = ModelClass(
                project_id=proj.id,
                subject_id=subj_db_id,
                sample_id=sample_db_id,
                category=category,
                subcategory=subcategory,
                platform=platform,
                modality=modality,
                measurement=measurement,
                resolution=resolution,
                spatial=spatial,
                uri=uri,
                format=fmt,
                status=status,
                file_role=file_role,
                validated=validated,
                n_obs=n_obs,
                n_vars=n_vars,
                phase=phase,
                size_mb=size_mb,
                md5=md5,
                **filtered_kwargs,
            )
            sess.add(obj)
            sess.commit()
            sess.refresh(obj)
            logger.info(
                "Registered BioData[%s] id=%d for project '%s' at %s",
                category,
                obj.id,
                project_name,
                uri,
            )
            return obj.id

    def get_biodata(self, biodata_id: int) -> dict[str, Any] | None:
        """Return a BioData dict (with type-specific columns) by id."""
        with self._session() as sess:
            row = sess.get(self._BioData, biodata_id)
            if row is None:
                return None
            # _to_dict via mapper.column_attrs gets all inherited columns
            return self._to_dict(row)

    def list_biodata(
        self,
        project_name: str | None = None,
        category: str | None = None,
        platform: str | None = None,
        sample_db_id: int | None = None,
        phase: str | None = None,
    ) -> list[dict[str, Any]]:
        """List BioData objects, optionally filtered."""
        with self._session() as sess:
            q = sess.query(self._BioData)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                q = q.filter_by(project_id=proj.id)
            if category is not None:
                q = q.filter_by(type=category)
            if platform is not None:
                q = q.filter_by(platform=platform)
            if sample_db_id is not None:
                q = q.filter_by(sample_id=sample_db_id)
            if phase is not None:
                q = q.filter_by(phase=phase)
            return [self._to_dict(r) for r in q.order_by(self._BioData.id).all()]

    def project_data_summary(self, project_name: str) -> dict[str, Any]:
        """Return counts of BioData objects grouped by category, modality, and platform."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return {
                    "project": project_name,
                    "total": 0,
                    "by_category": {},
                    "by_modality": {},
                    "by_platform": {},
                }
            rows = sess.query(self._BioData).filter_by(project_id=proj.id).all()
            by_cat: dict[str, int] = {}
            by_mod: dict[str, int] = {}
            by_plat: dict[str, int] = {}
            for r in rows:
                cat = r.type or "unknown"
                by_cat[cat] = by_cat.get(cat, 0) + 1
                mod = r.modality or "unknown"
                by_mod[mod] = by_mod.get(mod, 0) + 1
                plat = r.platform or "unknown"
                by_plat[plat] = by_plat.get(plat, 0) + 1
            return {
                "project": project_name,
                "total": len(rows),
                "by_category": by_cat,
                "by_modality": by_mod,
                "by_platform": by_plat,
            }

    # ------------------------------------------------------------------
    # Migration: datasets -> biodata
    # ------------------------------------------------------------------

    def migrate_datasets_to_biodata(self) -> int:
        """Migrate existing Dataset rows to BioData. Returns count migrated.

        Infers BioData subtype from project platform, file format, and file_role.
        Sets ``legacy_dataset_id`` on the new BioData row and ``bio_data_id`` on
        the Dataset row for bidirectional linking.
        """
        count = 0
        with self._session() as sess:
            datasets = sess.query(self._Dataset).all()
            for ds in datasets:
                # Skip if already migrated
                if getattr(ds, "bio_data_id", None) is not None:
                    continue
                proj = sess.get(self._Project, ds.project_id)
                proj_platform = proj.platform if proj else None

                # Infer subtype using platform registry when available
                image_type = None
                try:
                    from sc_tools.biodata import platform_for_project

                    pspec = platform_for_project(proj_platform) if proj_platform else None
                except ImportError:
                    pspec = None

                # Image files (tiff/image role) are always BioImage
                if ds.format == "tiff" or ds.file_role == "image":
                    ModelClass = self._BioImage
                    category = "image"
                    measurement = "protein" if proj_platform in ("imc", "mibi") else "morphology"
                    image_type = "multiplexed" if proj_platform in ("imc", "mibi") else None
                elif pspec is not None:
                    # Use platform registry to classify
                    type_map = {
                        "spatial_seq": self._SpatialSeqData,
                        "image": self._BioImage,
                        "rnaseq": self._RNASeqData,
                        "epigenomics": self._EpigenomicsData,
                        "genome_seq": self._GenomeSeqData,
                    }
                    ModelClass = type_map.get(pspec.biodata_type, self._SpatialSeqData)
                    category = pspec.biodata_type
                    measurement = pspec.measurement
                    if category == "image":
                        image_type = "multiplexed" if pspec.subcategory == "multiplexed" else None
                else:
                    # Unknown platform — default to SpatialSeqData (most existing data)
                    ModelClass = self._SpatialSeqData
                    category = "spatial_seq"
                    measurement = "rna"

                kwargs: dict[str, Any] = {
                    "project_id": ds.project_id,
                    "category": category,
                    "platform": proj_platform,
                    "measurement": measurement,
                    "uri": ds.uri,
                    "format": ds.format,
                    "status": ds.status,
                    "file_role": ds.file_role or "primary",
                    "validated": ds.validated or False,
                    "n_obs": ds.n_obs,
                    "n_vars": ds.n_vars,
                    "phase": ds.phase,
                    "size_mb": ds.size_mb,
                    "md5": ds.md5,
                    "legacy_dataset_id": ds.id,
                }
                if image_type and ModelClass == self._BioImage:
                    kwargs["image_type"] = image_type

                obj = ModelClass(**kwargs)
                sess.add(obj)
                sess.flush()
                # Forward-link
                ds.bio_data_id = obj.id
                count += 1

            sess.commit()
        logger.info("Migrated %d datasets to BioData", count)
        return count

    # ------------------------------------------------------------------
    # Status summary
    # ------------------------------------------------------------------

    def status(self) -> dict[str, Any]:
        """Return a high-level summary of current registry state."""
        with self._session() as sess:
            n_projects = sess.query(self._Project).count()
            n_datasets = sess.query(self._Dataset).count()
            n_subjects = sess.query(self._Subject).count()
            n_samples = sess.query(self._Sample).count()
            n_biodata = sess.query(self._BioData).count()
            active_jobs = (
                sess.query(self._SlurmJob)
                .filter(self._SlurmJob.status.in_(["submitted", "running"]))
                .count()
            )
            running_tasks = sess.query(self._AgentTask).filter_by(status="running").count()
            projects = [p.name for p in sess.query(self._Project).filter_by(status="active").all()]
            # Per-project phase summary
            phase_summary: dict[str, dict[str, int]] = {}
            for proj_name in projects:
                proj = sess.query(self._Project).filter_by(name=proj_name).first()
                if proj is None:
                    continue
                rows = sess.query(self._ProjectPhase).filter_by(project_id=proj.id).all()
                counts: dict[str, int] = {}
                for r in rows:
                    counts[r.status] = counts.get(r.status, 0) + 1
                if counts:
                    phase_summary[proj_name] = counts
        return {
            "n_projects": n_projects,
            "n_datasets": n_datasets,
            "n_subjects": n_subjects,
            "n_samples": n_samples,
            "n_biodata": n_biodata,
            "active_slurm_jobs": active_jobs,
            "running_agent_tasks": running_tasks,
            "active_projects": projects,
            "phase_summary": phase_summary,
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_dict(row: Any) -> dict[str, Any]:
        """Convert an ORM row to a dict.

        Uses ``mapper.column_attrs`` instead of ``__table__.columns`` so
        that Joined Table Inheritance (JTI) subclasses include all columns
        from both parent and child tables.
        """
        from sqlalchemy import inspect as sa_inspect

        mapper = sa_inspect(type(row))
        return {attr.key: getattr(row, attr.key) for attr in mapper.column_attrs}


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------


def _cli_status() -> None:
    """Print a registry status summary to stdout."""
    reg = Registry()
    s = reg.status()
    print("\nsc_tools Registry Status")
    print("=" * 40)
    print(f"  Projects          : {s['n_projects']}")
    print(f"  Datasets          : {s['n_datasets']}")
    print(f"  Subjects          : {s['n_subjects']}")
    print(f"  Samples           : {s['n_samples']}")
    print(f"  BioData objects   : {s['n_biodata']}")
    print(f"  Active SLURM jobs : {s['active_slurm_jobs']}")
    print(f"  Running tasks     : {s['running_agent_tasks']}")
    if s["active_projects"]:
        print("\n  Active projects:")
        for proj_name in s["active_projects"]:
            phases_done: list[str] = []
            proj = reg.get_project(proj_name)
            if proj:
                phases_done = json.loads(proj.get("phases_complete") or "[]")
            label = f" [phases: {', '.join(phases_done)}]" if phases_done else ""
            print(f"    - {proj_name}{label}")
            # Per-phase status
            phase_summary = s.get("phase_summary", {}).get(proj_name, {})
            if phase_summary:
                parts = ", ".join(f"{k}={v}" for k, v in sorted(phase_summary.items()))
                print(f"      phases: {parts}")
            datasets = reg.list_datasets(project_name=proj_name)
            for ds in datasets[-5:]:
                role = ds.get("file_role", "primary")
                print(
                    f"        phase={ds['phase']} role={role} status={ds['status']} uri={ds['uri']}"
                )
    print()


__all__ = ["Registry"]
