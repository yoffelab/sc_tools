"""Data registry for sc_tools projects.

Tracks data objects (AnnData checkpoints, images, etc.) and patients across
all active projects.  Uses **SQLite** by default (zero-config) with optional
**PostgreSQL** support via the ``SC_TOOLS_REGISTRY_URL`` environment variable.

Quick start
-----------
::

    from sc_tools.registry import Registry

    reg = Registry()
    reg.add_project("ggo_visium", platform="visium",
                    domain="spatial_transcriptomics")
    ds_id = reg.register_data(
        "ggo_visium", phase="qc_filter",
        uri="sftp://brb//athena/.../adata.raw.h5ad",
        fmt="h5ad",
        file_role="primary",
    )

CLI
---
::

    python -m sc_tools registry status

Environment
-----------
``SC_TOOLS_REGISTRY_URL``
    Optional SQLAlchemy-compatible DB URL.  Defaults to
    ``sqlite:///~/.sc_tools/registry.db``.

Phase slugs
-----------
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


_URI_SCHEMES = ("sftp://", "s3://", "gs://", "az://", "http://", "https://", "phase://")


def _validate_uri(uri: str) -> str:
    """Ensure URI is an absolute path or has a recognized scheme."""
    if any(uri.startswith(s) for s in _URI_SCHEMES):
        return uri
    if os.path.isabs(uri):
        return uri
    raise ValueError(
        f"URI must be an absolute path or use a scheme ({', '.join(_URI_SCHEMES)}). Got: {uri!r}"
    )


# ---------------------------------------------------------------------------
# ORM models — built lazily to avoid top-level SQLAlchemy import errors
# ---------------------------------------------------------------------------


def _build_models(Base: Any) -> tuple:
    """Return ORM model classes.

    Returns a 6-tuple:
        Project, Data, Patient, DataSource, ProjectDataSource, PatientDataMap
    """
    from sqlalchemy import (
        Column,
        Float,
        ForeignKey,
        Integer,
        String,
        Text,
    )
    from sqlalchemy.orm import mapped_column, relationship

    # Use JSONB on PostgreSQL, JSON on SQLite
    try:
        from sqlalchemy.dialects.postgresql import JSONB as JSONType
    except ImportError:
        from sqlalchemy import JSON as JSONType

    class Project(Base):
        __tablename__ = "projects"

        id = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String, unique=True, nullable=False)
        platform = Column(String)
        domain = Column(String)
        status = Column(String, default="active")
        created_at = Column(String, default=_utcnow)

        data = relationship("Data", back_populates="project", cascade="all, delete-orphan")
        data_source_links = relationship(
            "ProjectDataSource", back_populates="project", cascade="all, delete-orphan"
        )

    class Data(Base):
        __tablename__ = "data_processing_phase"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(
            Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False
        )
        phase = Column(String)
        status = Column(String, default="pending")
        uri = Column(Text, nullable=False)
        format = Column(String)
        platform = Column(String)
        category = Column(String)
        file_role = Column(String, default="primary")
        n_obs = Column(Integer)
        n_vars = Column(Integer)
        size_mb = Column(Float)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        project = relationship("Project", back_populates="data")
        patient_links = relationship("PatientDataMap", back_populates="data", cascade="all, delete-orphan")

    class Patient(Base):
        __tablename__ = "patients"

        id = Column(Integer, primary_key=True, autoincrement=True)
        patient_id = Column(String, unique=True, nullable=False)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)

        data_links = relationship("PatientDataMap", back_populates="patient", cascade="all, delete-orphan")

    class DataSource(Base):
        """A raw data source -- HPC directory, public dataset, or external repository."""

        __tablename__ = "data_inventory"

        id = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String, unique=True, nullable=False)
        description = Column(Text)
        uri = Column(Text, nullable=False)
        platform = Column(String)
        domain = Column(String)
        imaging_modality = Column(String)
        source_type = Column(String)
        organism = Column(String)
        tissue = Column(String)
        disease = Column(String)
        n_samples = Column(Integer)
        n_cells = Column(Integer)
        publication = Column(String)
        access_notes = Column(Text)
        status = Column(String, default="available")
        created_at = Column(String, default=_utcnow)

        project_links = relationship(
            "ProjectDataSource", back_populates="data_source", cascade="all, delete-orphan"
        )

    class ProjectDataSource(Base):
        """Many-to-many join between Projects and DataSources."""

        __tablename__ = "data_project_map"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(
            Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False
        )
        data_source_id = Column(
            Integer, ForeignKey("data_inventory.id", ondelete="CASCADE"), nullable=False
        )
        role = Column(String, default="input")
        notes = Column(Text)

        project = relationship("Project", back_populates="data_source_links")
        data_source = relationship("DataSource", back_populates="project_links")

    class PatientDataMap(Base):
        """Many-to-many join between Patients and Data."""

        __tablename__ = "patient_data_map"

        id = Column(Integer, primary_key=True, autoincrement=True)
        patient_id = Column(
            Integer, ForeignKey("patients.id", ondelete="CASCADE"), nullable=False
        )
        data_id = Column(
            Integer, ForeignKey("data_processing_phase.id", ondelete="CASCADE"), nullable=False
        )
        role = Column(String, default="source")
        notes = Column(Text)

        patient = relationship("Patient", back_populates="data_links")
        data = relationship("Data", back_populates="patient_links")

    return (
        Project,
        Data,
        Patient,
        DataSource,
        ProjectDataSource,
        PatientDataMap,
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
            self._Data,
            self._Patient,
            self._DataSource,
            self._ProjectDataSource,
            self._PatientDataMap,
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
        data_type: str | None = None,  # kept for backward compat, ignored
        *,
        domain: str | None = None,
        imaging_modality: str | None = None,  # ignored (column dropped)
        project_type: str = "internal",  # ignored (column dropped)
        visibility: str = "private",  # ignored (column dropped)
    ) -> int:
        """Register a new project. Returns project id.

        If a project with the same name exists its id is returned without
        creating a duplicate.
        """
        with self._session() as sess:
            existing = sess.query(self._Project).filter_by(name=name).first()
            if existing:
                logger.info("Project '%s' already registered (id=%d)", name, existing.id)
                return existing.id
            proj = self._Project(
                name=name,
                platform=platform,
                domain=domain,
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

    def mark_phase_complete(self, project_name: str, phase: str) -> None:
        """Mark a pipeline phase as complete.

        Updates all ``data`` rows for this project+phase to status='ready'.
        """
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")
            rows = (
                sess.query(self._Data)
                .filter_by(project_id=proj.id, phase=phase)
                .all()
            )
            now = _utcnow()
            for r in rows:
                r.status = "ready"
                r.updated_at = now
            sess.commit()

    # ------------------------------------------------------------------
    # Data (replaces datasets + bio_data)
    # ------------------------------------------------------------------

    def register_data(
        self,
        project_name: str,
        phase: str,
        uri: str,
        *,
        fmt: str | None = "h5ad",
        platform: str | None = None,
        category: str | None = None,
        status: str = "ready",
        file_role: str = "primary",
        n_obs: int | None = None,
        n_vars: int | None = None,
        size_mb: float | None = None,
        patient_id: str | None = None,
        metadata: dict | None = None,
        sample_id: str | None = None,  # backward compat, ignored
    ) -> int:
        """Register a data object. Returns data id.

        Parameters
        ----------
        project_name:
            Project name (must already exist).
        phase:
            Phase slug (e.g. ``"qc_filter"``, ``"scoring"``).
        uri:
            Path or URI to the data file.
        fmt:
            File format: ``"h5ad"`` (default), ``"zarr"``, ``"tiff"``, ``"tsv"``.
        platform:
            Platform slug. If None, inferred from project.
        category:
            Data category. If None, inferred from platform.
        file_role:
            Role: ``"primary"`` | ``"supplementary"`` | ``"entry_point"`` etc.
        patient_id:
            De-identified patient ID string (optional). Must already exist.
        metadata:
            Arbitrary JSONB metadata dict.
        """
        _validate_uri(uri)
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found. Call add_project() first.")

            # Infer platform/category from project if not provided
            if platform is None:
                platform = proj.platform
            if category is None:
                category = "spatial_seq"  # default
                try:
                    from sc_tools.biodata import platform_for_project

                    pspec = platform_for_project(platform) if platform else None
                    if pspec:
                        category = pspec.biodata_type
                except (ImportError, KeyError):
                    pass

            # Auto-fill format from URI extension
            if fmt is None and uri:
                ext = uri.rsplit(".", 1)[-1].lower() if "." in uri else ""
                if ext in ("h5ad", "zarr", "tiff", "tsv", "fastq", "bam", "bed", "csv"):
                    fmt = ext

            now = _utcnow()
            obj = self._Data(
                project_id=proj.id,
                phase=phase,
                status=status,
                uri=uri,
                format=fmt,
                platform=platform,
                category=category,
                file_role=file_role,
                n_obs=n_obs,
                n_vars=n_vars,
                size_mb=size_mb,
                meta=metadata or {},
                created_at=now,
                updated_at=now,
            )
            sess.add(obj)
            sess.flush()

            # Link patient via join table if provided
            if patient_id is not None:
                pat = sess.query(self._Patient).filter_by(patient_id=patient_id).first()
                if pat is not None:
                    link = self._PatientDataMap(patient_id=pat.id, data_id=obj.id)
                    sess.add(link)

            sess.commit()
            sess.refresh(obj)
            logger.info(
                "Registered data id=%d for project '%s' phase=%s at %s",
                obj.id,
                project_name,
                phase,
                uri,
            )
            return obj.id

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
        """Backward-compatible alias for :meth:`register_data`."""
        meta: dict[str, Any] = {}
        if md5:
            meta["md5"] = md5
        if validated:
            meta["validated"] = True
        return self.register_data(
            project_name=project_name,
            phase=phase,
            uri=uri,
            fmt=fmt,
            status=status,
            file_role=file_role,
            n_obs=n_obs,
            n_vars=n_vars,
            size_mb=size_mb,
            metadata=meta if meta else None,
        )

    def register_biodata(
        self,
        project_name: str,
        category: str,
        platform: str,
        uri: str,
        *,
        fmt: str | None = None,
        status: str = "pending",
        file_role: str = "primary",
        phase: str | None = None,
        n_obs: int | None = None,
        n_vars: int | None = None,
        size_mb: float | None = None,
        md5: str | None = None,
        subject_id: str | None = None,
        sample_db_id: int | None = None,
        **type_kwargs: Any,
    ) -> int:
        """Backward-compatible alias for :meth:`register_data`.

        Accepts the old BioData-style arguments and routes to register_data.
        """
        meta: dict[str, Any] = {}
        if md5:
            meta["md5"] = md5
        if type_kwargs:
            meta.update(type_kwargs)
        return self.register_data(
            project_name=project_name,
            phase=phase or "unknown",
            uri=uri,
            fmt=fmt,
            platform=platform,
            category=category,
            status=status,
            file_role=file_role,
            n_obs=n_obs,
            n_vars=n_vars,
            size_mb=size_mb,
            patient_id=subject_id,
            metadata=meta if meta else None,
        )

    def get_dataset_uri(
        self,
        project_name: str,
        phase: str,
        sample_id: str | None = None,
    ) -> str | None:
        """Return the URI for a checkpoint (most recently registered)."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return None
            q = sess.query(self._Data).filter_by(project_id=proj.id, phase=phase)
            row = q.order_by(self._Data.id.desc()).first()
            return row.uri if row else None

    def list_datasets(
        self,
        project_name: str | None = None,
        phase: str | None = None,
    ) -> list[dict[str, Any]]:
        """List data entries, optionally filtered by project and/or phase."""
        with self._session() as sess:
            q = sess.query(self._Data)
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
        """Update data row status."""
        with self._session() as sess:
            row = sess.get(self._Data, dataset_id)
            if row is None:
                raise ValueError(f"Data id={dataset_id} not found")
            row.status = status
            row.updated_at = _utcnow()
            sess.commit()

    # ------------------------------------------------------------------
    # Phase queries (derived from data table)
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
        """Create or update a data row for phase tracking.

        For backward compatibility with the old project_phases table.
        Creates a metadata-only data row if no data row exists for this
        project+phase, or updates the metadata of the most recent one.
        """
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")

            existing = (
                sess.query(self._Data)
                .filter_by(project_id=proj.id, phase=phase)
                .order_by(self._Data.id.desc())
                .first()
            )
            now = _utcnow()
            if existing is not None:
                existing.status = status
                existing.updated_at = now
                meta = dict(existing.meta or {})
                if n_obs is not None:
                    existing.n_obs = n_obs
                if n_vars is not None:
                    existing.n_vars = n_vars
                if n_samples is not None:
                    meta["n_samples"] = n_samples
                if notes is not None:
                    meta["notes"] = notes
                if entry_phase:
                    meta["entry_phase"] = True
                if primary_dataset_id is not None:
                    meta["primary_dataset_id"] = primary_dataset_id
                existing.meta = meta
            else:
                meta: dict[str, Any] = {}
                if n_samples is not None:
                    meta["n_samples"] = n_samples
                if notes is not None:
                    meta["notes"] = notes
                if entry_phase:
                    meta["entry_phase"] = True
                if primary_dataset_id is not None:
                    meta["primary_dataset_id"] = primary_dataset_id
                row = self._Data(
                    project_id=proj.id,
                    phase=phase,
                    status=status,
                    uri=f"phase://{project_name}/{phase}",
                    file_role="phase_marker",
                    n_obs=n_obs,
                    n_vars=n_vars,
                    platform=proj.platform,
                    meta=meta,
                    created_at=now,
                    updated_at=now,
                )
                sess.add(row)
            sess.commit()

    def get_phase(self, project_name: str, phase: str) -> dict[str, Any] | None:
        """Return phase info derived from data table."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return None
            row = (
                sess.query(self._Data)
                .filter_by(project_id=proj.id, phase=phase)
                .order_by(self._Data.id.desc())
                .first()
            )
            if row is None:
                return None
            meta = row.meta or {}
            return {
                "phase": phase,
                "status": row.status,
                "entry_phase": meta.get("entry_phase", False),
                "n_obs": row.n_obs,
                "n_vars": row.n_vars,
                "n_samples": meta.get("n_samples"),
                "primary_dataset_id": meta.get("primary_dataset_id"),
                "notes": meta.get("notes"),
                "updated_at": row.updated_at,
            }

    def list_phases(self, project_name: str) -> list[dict[str, Any]]:
        """Return phase summaries for a project, derived from data table."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return []
            from sqlalchemy import func

            # Get distinct phases with their latest status
            subq = (
                sess.query(
                    self._Data.phase,
                    func.max(self._Data.id).label("max_id"),
                )
                .filter_by(project_id=proj.id)
                .filter(self._Data.phase.isnot(None))
                .group_by(self._Data.phase)
                .subquery()
            )
            rows = (
                sess.query(self._Data)
                .join(subq, self._Data.id == subq.c.max_id)
                .order_by(self._Data.phase)
                .all()
            )
            result = []
            for r in rows:
                rmeta = r.meta or {}
                result.append(
                    {
                        "phase": r.phase,
                        "status": r.status,
                        "n_obs": r.n_obs,
                        "n_vars": r.n_vars,
                        "n_samples": rmeta.get("n_samples"),
                        "updated_at": r.updated_at,
                    }
                )
            return result

    # ------------------------------------------------------------------
    # Data sources (UNCHANGED)
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
        """Register a data source. Returns data source id. Idempotent."""
        _validate_uri(uri)
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
        """Link a project to a data source. Returns link id. Idempotent."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found.")
            src = sess.query(self._DataSource).filter_by(name=data_source_name).first()
            if src is None:
                raise ValueError(f"DataSource '{data_source_name}' not found.")
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
                "Linked project '%s' -> data source '%s' (role=%s)",
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

    # ------------------------------------------------------------------
    # Patients (replaces subjects + samples)
    # ------------------------------------------------------------------

    def add_patient(
        self,
        patient_id: str,
        metadata: dict | None = None,
        **clinical_kwargs: Any,
    ) -> int:
        """Register a patient. Returns patient DB id.

        Idempotent: if a patient with the same ``patient_id`` exists, its
        DB id is returned without creating a duplicate.

        Parameters
        ----------
        patient_id:
            Unique de-identified identifier (e.g. ``"PT001"``).
        metadata:
            JSONB metadata dict with clinical info (organism, sex, diagnosis, etc.).
        **clinical_kwargs:
            Individual clinical fields merged into metadata dict.
        """
        combined_meta = dict(metadata or {})
        combined_meta.update(clinical_kwargs)

        with self._session() as sess:
            existing = sess.query(self._Patient).filter_by(patient_id=patient_id).first()
            if existing:
                logger.info("Patient '%s' already registered (id=%d)", patient_id, existing.id)
                return existing.id
            pat = self._Patient(
                patient_id=patient_id,
                meta=combined_meta,
            )
            sess.add(pat)
            sess.commit()
            sess.refresh(pat)
            logger.info("Registered patient '%s' (id=%d)", patient_id, pat.id)
            return pat.id

    def add_subject(
        self,
        subject_id: str,
        organism: str = "human",
        **clinical_kwargs: Any,
    ) -> int:
        """Backward-compatible alias for :meth:`add_patient`.

        Maps subject_id -> patient_id and packs clinical columns into metadata.
        """
        meta = {"organism": organism}
        meta.update(clinical_kwargs)
        return self.add_patient(subject_id, metadata=meta)

    def get_patient(self, patient_id: str) -> dict[str, Any] | None:
        """Return patient dict by de-identified ID, or None."""
        with self._session() as sess:
            row = sess.query(self._Patient).filter_by(patient_id=patient_id).first()
            return self._to_dict(row) if row else None

    # Backward compat alias
    def get_subject(self, subject_id: str) -> dict[str, Any] | None:
        """Backward-compatible alias for :meth:`get_patient`."""
        return self.get_patient(subject_id)

    def list_patients(
        self,
        diagnosis: str | None = None,
        tissue: str | None = None,
    ) -> list[dict[str, Any]]:
        """Return patients, optionally filtered by metadata fields."""
        with self._session() as sess:
            q = sess.query(self._Patient)
            # Filter using JSONB operators on PostgreSQL
            if diagnosis is not None:
                q = q.filter(
                    self._Patient.meta["diagnosis"].astext.ilike(f"%{diagnosis}%")
                )
            if tissue is not None:
                q = q.filter(
                    self._Patient.meta["tissue_of_origin"].astext.ilike(f"%{tissue}%")
                )
            return [self._to_dict(r) for r in q.order_by(self._Patient.patient_id).all()]

    def list_subjects(
        self,
        project_name: str | None = None,
        diagnosis: str | None = None,
        tissue: str | None = None,
    ) -> list[dict[str, Any]]:
        """Backward-compatible alias for :meth:`list_patients`.

        The project_name filter now checks for patients that have data in that project.
        """
        with self._session() as sess:
            q = sess.query(self._Patient)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                # Find patients linked to data in this project via join table
                patient_ids_subq = (
                    sess.query(self._PatientDataMap.patient_id)
                    .join(self._Data, self._PatientDataMap.data_id == self._Data.id)
                    .filter(self._Data.project_id == proj.id)
                    .distinct()
                    .scalar_subquery()
                )
                q = q.filter(self._Patient.id.in_(patient_ids_subq))
            if diagnosis is not None:
                q = q.filter(
                    self._Patient.meta["diagnosis"].astext.ilike(f"%{diagnosis}%")
                )
            if tissue is not None:
                q = q.filter(
                    self._Patient.meta["tissue_of_origin"].astext.ilike(f"%{tissue}%")
                )
            results = []
            for r in q.order_by(self._Patient.patient_id).all():
                d = self._to_dict(r)
                # Flatten metadata for backward compat
                meta = d.get("meta") or {}
                d["subject_id"] = d["patient_id"]
                d["organism"] = meta.get("organism")
                d["sex"] = meta.get("sex")
                d["diagnosis"] = meta.get("diagnosis")
                d["tissue_of_origin"] = meta.get("tissue_of_origin")
                results.append(d)
            return results

    # ------------------------------------------------------------------
    # BioData queries (backward compat)
    # ------------------------------------------------------------------

    def list_biodata(
        self,
        project_name: str | None = None,
        category: str | None = None,
        platform: str | None = None,
        sample_db_id: int | None = None,
        phase: str | None = None,
    ) -> list[dict[str, Any]]:
        """List data objects. Backward-compatible alias for list_datasets."""
        with self._session() as sess:
            q = sess.query(self._Data)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                q = q.filter_by(project_id=proj.id)
            if category is not None:
                q = q.filter_by(category=category)
            if platform is not None:
                q = q.filter_by(platform=platform)
            if phase is not None:
                q = q.filter_by(phase=phase)
            results = []
            for r in q.order_by(self._Data.id).all():
                d = self._to_dict(r)
                # Provide backward-compat 'type' field from meta or category
                meta = d.get("meta") or {}
                d["type"] = meta.get("type", d.get("category", "unknown"))
                d["modality"] = meta.get("modality")
                results.append(d)
            return results

    def get_biodata(self, biodata_id: int) -> dict[str, Any] | None:
        """Return a data dict by id."""
        with self._session() as sess:
            row = sess.get(self._Data, biodata_id)
            if row is None:
                return None
            return self._to_dict(row)

    def project_data_summary(self, project_name: str) -> dict[str, Any]:
        """Return counts of data objects grouped by category and platform."""
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
            rows = sess.query(self._Data).filter_by(project_id=proj.id).all()
            by_cat: dict[str, int] = {}
            by_mod: dict[str, int] = {}
            by_plat: dict[str, int] = {}
            for r in rows:
                cat = r.category or "unknown"
                by_cat[cat] = by_cat.get(cat, 0) + 1
                rmeta = r.meta or {}
                mod = rmeta.get("modality", "unknown")
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
    # Status summary
    # ------------------------------------------------------------------

    def status(self) -> dict[str, Any]:
        """Return a high-level summary of current registry state."""
        with self._session() as sess:
            n_projects = sess.query(self._Project).count()
            n_data = sess.query(self._Data).count()
            n_patients = sess.query(self._Patient).count()
            projects = [p.name for p in sess.query(self._Project).filter_by(status="active").all()]

            # Per-project phase summary
            phase_summary: dict[str, dict[str, int]] = {}
            for proj_name in projects:
                proj = sess.query(self._Project).filter_by(name=proj_name).first()
                if proj is None:
                    continue
                from sqlalchemy import func

                phase_counts = (
                    sess.query(self._Data.status, func.count())
                    .filter_by(project_id=proj.id)
                    .filter(self._Data.phase.isnot(None))
                    .group_by(self._Data.status)
                    .all()
                )
                counts = {status: count for status, count in phase_counts}
                if counts:
                    phase_summary[proj_name] = counts

        return {
            "n_projects": n_projects,
            "n_datasets": n_data,
            "n_data": n_data,
            "n_patients": n_patients,
            # Backward-compat keys
            "n_subjects": n_patients,
            "n_samples": 0,
            "n_biodata": n_data,
            "active_slurm_jobs": 0,
            "running_agent_tasks": 0,
            "active_projects": projects,
            "phase_summary": phase_summary,
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_dict(row: Any) -> dict[str, Any]:
        """Convert an ORM row to a dict."""
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
    print(f"  Data objects      : {s['n_data']}")
    print(f"  Patients          : {s['n_patients']}")
    if s["active_projects"]:
        print("\n  Active projects:")
        for proj_name in s["active_projects"]:
            print(f"    - {proj_name}")
            phase_summary = s.get("phase_summary", {}).get(proj_name, {})
            if phase_summary:
                parts = ", ".join(f"{k}={v}" for k, v in sorted(phase_summary.items()))
                print(f"      data status: {parts}")
            datasets = reg.list_datasets(project_name=proj_name)
            for ds in datasets[-5:]:
                role = ds.get("file_role", "primary")
                print(
                    f"        phase={ds['phase']} role={role} status={ds['status']} uri={ds['uri']}"
                )
    print()


__all__ = ["Registry"]
