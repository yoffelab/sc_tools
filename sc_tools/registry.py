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


def _build_models(Base: Any) -> tuple[Any, Any, Any, Any, Any, Any, Any]:
    """Return ORM model classes (Project, Dataset, SlurmJob, AgentTask, ProjectPhase,
    DataSource, ProjectDataSource)."""
    from sqlalchemy import Boolean, Column, Float, ForeignKey, Integer, String, Text
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

    return Project, Dataset, SlurmJob, AgentTask, ProjectPhase, DataSource, ProjectDataSource


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
            logger.info("Registered dataset %s/%s at %s (id=%d)", project_name, phase, uri, ds.id)
            return ds.id

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
    # Status summary
    # ------------------------------------------------------------------

    def status(self) -> dict[str, Any]:
        """Return a high-level summary of current registry state."""
        with self._session() as sess:
            n_projects = sess.query(self._Project).count()
            n_datasets = sess.query(self._Dataset).count()
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
        return {col.name: getattr(row, col.name) for col in row.__table__.columns}


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
