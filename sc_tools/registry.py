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

import logging
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
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
# ORM models -- built lazily to avoid top-level SQLAlchemy import errors
# ---------------------------------------------------------------------------


def _build_models(Base: Any) -> SimpleNamespace:
    """Return ORM model classes as a SimpleNamespace.

    Attributes
    ----------
    Project, DataSource, InventoryItem, Dataset, DatasetMember,
    ProjectDataset, ProjectPhase, Patient, Sample, Provenance
    """
    from sqlalchemy import (
        JSON as JSONType,
    )
    from sqlalchemy import (
        Boolean,
        Column,
        Float,
        ForeignKey,
        Integer,
        String,
        Text,
    )
    from sqlalchemy.orm import relationship

    # ------------------------------------------------------------------
    # Existing tables (unchanged)
    # ------------------------------------------------------------------

    class Project(Base):
        __tablename__ = "projects"

        id = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String, unique=True, nullable=False)
        status = Column(String, default="active")
        created_at = Column(String, default=_utcnow)

        phase_links = relationship(
            "ProjectPhase", back_populates="project", cascade="all, delete-orphan"
        )
        dataset_links = relationship(
            "ProjectDataset", back_populates="project", cascade="all, delete-orphan"
        )

    class Patient(Base):
        __tablename__ = "patients"

        id = Column(Integer, primary_key=True, autoincrement=True)
        patient_id = Column(String, unique=True, nullable=False)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)

        sample_links = relationship(
            "Sample", back_populates="patient", cascade="all, delete-orphan"
        )

    class Sample(Base):
        """A biological specimen collected from a patient."""

        __tablename__ = "samples"

        id = Column(Integer, primary_key=True, autoincrement=True)
        patient_id = Column(
            Integer, ForeignKey("patients.id", ondelete="SET NULL"), nullable=True
        )
        sample_id = Column(String, unique=True, nullable=False)
        tissue = Column(String)
        collection_date = Column(String)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)

        patient = relationship("Patient", back_populates="sample_links")
        inventory_items = relationship("InventoryItem", back_populates="sample")

    # ------------------------------------------------------------------
    # Layer 0: data_sources
    # ------------------------------------------------------------------

    class DataSource(Base):
        """A raw data source -- HPC directory, public dataset, or external repository."""

        __tablename__ = "data_sources"

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
        status = Column(String, default="discovered")
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)

        inventory_items = relationship("InventoryItem", back_populates="data_source")

    # ------------------------------------------------------------------
    # New Layer 1: inventory_items
    # ------------------------------------------------------------------

    class InventoryItem(Base):
        """A clean, ingested data file (h5ad) produced from a data source."""

        __tablename__ = "inventory_items"

        id = Column(Integer, primary_key=True, autoincrement=True)
        data_source_id = Column(
            Integer, ForeignKey("data_sources.id", ondelete="SET NULL"), nullable=True
        )
        sample_id = Column(
            Integer, ForeignKey("samples.id", ondelete="SET NULL"), nullable=True
        )
        name = Column(String, unique=True, nullable=False)
        uri = Column(Text, nullable=False)
        modality = Column(String, nullable=False)
        platform = Column(String)
        format = Column(String, default="h5ad")
        n_obs = Column(Integer)
        n_vars = Column(Integer)
        size_mb = Column(Float)
        organism = Column(String)
        tissue = Column(String)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        data_source = relationship("DataSource", back_populates="inventory_items")
        sample = relationship("Sample", back_populates="inventory_items")
        dataset_memberships = relationship("DatasetMember", back_populates="inventory_item")

    # ------------------------------------------------------------------
    # New Layer 2: datasets + dataset_members
    # ------------------------------------------------------------------

    class Dataset(Base):
        """A versioned assembly of inventory items (MuData or AnnData)."""

        __tablename__ = "datasets"

        id = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String, nullable=False)
        version = Column(Integer, nullable=False, default=1)
        description = Column(Text)
        uri = Column(Text)
        format = Column(String, default="mudata")
        n_obs = Column(Integer)
        size_mb = Column(Float)
        is_current = Column(Boolean, default=True)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        members = relationship(
            "DatasetMember", back_populates="dataset", cascade="all, delete-orphan"
        )
        project_links = relationship(
            "ProjectDataset", back_populates="dataset", cascade="all, delete-orphan"
        )
        phases = relationship("ProjectPhase", back_populates="dataset")

    class DatasetMember(Base):
        """Composition link: which inventory items belong to a dataset."""

        __tablename__ = "dataset_members"

        id = Column(Integer, primary_key=True, autoincrement=True)
        dataset_id = Column(Integer, ForeignKey("datasets.id", ondelete="CASCADE"), nullable=False)
        inventory_id = Column(
            Integer, ForeignKey("inventory_items.id", ondelete="RESTRICT"), nullable=False
        )
        modality_key = Column(String, nullable=False)

        dataset = relationship("Dataset", back_populates="members")
        inventory_item = relationship("InventoryItem", back_populates="dataset_memberships")

    # ------------------------------------------------------------------
    # New Layer 3: project_datasets + project_phases
    # ------------------------------------------------------------------

    class ProjectDataset(Base):
        """Many-to-many join between Projects and Datasets."""

        __tablename__ = "project_datasets"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
        dataset_id = Column(Integer, ForeignKey("datasets.id", ondelete="CASCADE"), nullable=False)
        role = Column(String, default="primary")
        notes = Column(Text)

        project = relationship("Project", back_populates="dataset_links")
        dataset = relationship("Dataset", back_populates="project_links")

    class ProjectPhase(Base):
        """Phase tracking per project + dataset combination."""

        __tablename__ = "project_phases"

        id = Column(Integer, primary_key=True, autoincrement=True)
        project_id = Column(Integer, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
        dataset_id = Column(Integer, ForeignKey("datasets.id", ondelete="RESTRICT"), nullable=False)
        phase_group = Column(String, nullable=False)
        subphase = Column(String, nullable=False)
        status = Column(String, default="pending")
        uri = Column(Text)
        n_obs = Column(Integer)
        n_vars = Column(Integer)
        meta = Column("metadata", JSONType, default=dict)
        created_at = Column(String, default=_utcnow)
        updated_at = Column(String, default=_utcnow)

        project = relationship("Project", back_populates="phase_links")
        dataset = relationship("Dataset", back_populates="phases")

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    class Provenance(Base):
        """Tracks tool versions, params, and environment for transformations."""

        __tablename__ = "provenance"

        id = Column(Integer, primary_key=True, autoincrement=True)
        target_type = Column(String, nullable=False)
        phase_id = Column(
            Integer, ForeignKey("project_phases.id", ondelete="SET NULL"), nullable=True
        )
        dataset_id = Column(Integer, ForeignKey("datasets.id", ondelete="SET NULL"), nullable=True)
        tool = Column(String, nullable=False)
        tool_version = Column(String)
        reference_genome = Column(String)
        reference_dataset = Column(Text)
        signature_source = Column(Text)
        n_input_obs = Column(Integer)
        n_output_obs = Column(Integer)
        params = Column(JSONType, default=dict)
        environment = Column(JSONType, default=dict)
        script_uri = Column(Text)
        agent = Column(String)
        created_at = Column(String, default=_utcnow)

        phase = relationship("ProjectPhase")
        dataset = relationship("Dataset")

    return SimpleNamespace(
        Project=Project,
        DataSource=DataSource,
        InventoryItem=InventoryItem,
        Dataset=Dataset,
        DatasetMember=DatasetMember,
        ProjectDataset=ProjectDataset,
        ProjectPhase=ProjectPhase,
        Patient=Patient,
        Sample=Sample,
        Provenance=Provenance,
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
        models = _build_models(Base)
        self._Project = models.Project
        self._Patient = models.Patient
        self._Sample = models.Sample
        self._DataSource = models.DataSource
        self._InventoryItem = models.InventoryItem
        self._Dataset = models.Dataset
        self._DatasetMember = models.DatasetMember
        self._ProjectDataset = models.ProjectDataset
        self._ProjectPhase = models.ProjectPhase
        self._Provenance = models.Provenance

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
        platform: str | None = None,  # ignored (column dropped)
        data_type: str | None = None,  # ignored
        *,
        domain: str | None = None,  # ignored (column dropped)
        imaging_modality: str | None = None,  # ignored
        project_type: str = "internal",  # ignored
        visibility: str = "private",  # ignored
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
            proj = self._Project(name=name)
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

    # ------------------------------------------------------------------
    # Inventory Items (NEW -- Layer 1)
    # ------------------------------------------------------------------

    def register_inventory_item(
        self,
        name: str,
        uri: str,
        modality: str,
        *,
        data_source_name: str | None = None,
        sample_name: str | None = None,
        platform: str | None = None,
        fmt: str = "h5ad",
        n_obs: int | None = None,
        n_vars: int | None = None,
        size_mb: float | None = None,
        organism: str | None = None,
        tissue: str | None = None,
        metadata: dict | None = None,
    ) -> int:
        """Register an inventory item. Returns item id. Idempotent by name."""
        _validate_uri(uri)
        with self._session() as sess:
            existing = sess.query(self._InventoryItem).filter_by(name=name).first()
            if existing:
                logger.info("InventoryItem '%s' already registered (id=%d)", name, existing.id)
                return existing.id

            data_source_id = None
            if data_source_name is not None:
                ds = sess.query(self._DataSource).filter_by(name=data_source_name).first()
                if ds is None:
                    raise ValueError(f"DataSource '{data_source_name}' not found.")
                data_source_id = ds.id

            sample_id = None
            if sample_name is not None:
                sample = sess.query(self._Sample).filter_by(sample_id=sample_name).first()
                if sample is None:
                    raise ValueError(f"Sample '{sample_name}' not found.")
                sample_id = sample.id

            now = _utcnow()
            item = self._InventoryItem(
                name=name,
                uri=uri,
                modality=modality,
                data_source_id=data_source_id,
                sample_id=sample_id,
                platform=platform,
                format=fmt,
                n_obs=n_obs,
                n_vars=n_vars,
                size_mb=size_mb,
                organism=organism,
                tissue=tissue,
                meta=metadata or {},
                created_at=now,
                updated_at=now,
            )
            sess.add(item)
            sess.commit()
            sess.refresh(item)
            logger.info("Registered inventory item '%s' (id=%d)", name, item.id)
            return item.id

    def get_inventory_item(self, name: str) -> dict[str, Any] | None:
        """Return inventory item dict by name, or None."""
        with self._session() as sess:
            row = sess.query(self._InventoryItem).filter_by(name=name).first()
            return self._to_dict(row) if row else None

    def list_inventory_items(
        self,
        modality: str | None = None,
        platform: str | None = None,
    ) -> list[dict[str, Any]]:
        """List inventory items with optional filters."""
        with self._session() as sess:
            q = sess.query(self._InventoryItem)
            if modality is not None:
                q = q.filter_by(modality=modality)
            if platform is not None:
                q = q.filter_by(platform=platform)
            return [self._to_dict(r) for r in q.order_by(self._InventoryItem.name).all()]

    # ------------------------------------------------------------------
    # Datasets (NEW -- Layer 2)
    # ------------------------------------------------------------------

    def create_dataset(
        self,
        name: str,
        *,
        description: str | None = None,
        fmt: str = "mudata",
    ) -> int:
        """Create a dataset. Returns dataset id. Idempotent for name+version=1."""
        with self._session() as sess:
            existing = sess.query(self._Dataset).filter_by(name=name, version=1).first()
            if existing:
                logger.info("Dataset '%s' v1 already exists (id=%d)", name, existing.id)
                return existing.id
            now = _utcnow()
            ds = self._Dataset(
                name=name,
                version=1,
                description=description,
                format=fmt,
                is_current=True,
                created_at=now,
                updated_at=now,
            )
            sess.add(ds)
            sess.commit()
            sess.refresh(ds)
            logger.info("Created dataset '%s' v1 (id=%d)", name, ds.id)
            return ds.id

    def get_dataset(
        self,
        name: str,
        version: int | None = None,
    ) -> dict[str, Any] | None:
        """Return dataset dict. If version is None, get is_current or highest version."""
        with self._session() as sess:
            if version is not None:
                row = sess.query(self._Dataset).filter_by(name=name, version=version).first()
            else:
                # Prefer is_current=True, fall back to highest version
                row = sess.query(self._Dataset).filter_by(name=name, is_current=True).first()
                if row is None:
                    row = (
                        sess.query(self._Dataset)
                        .filter_by(name=name)
                        .order_by(self._Dataset.version.desc())
                        .first()
                    )
            return self._to_dict(row) if row else None

    def add_dataset_member(
        self,
        dataset_name: str,
        inventory_name: str,
        modality_key: str,
    ) -> int:
        """Add an inventory item to a dataset. Returns member id."""
        with self._session() as sess:
            ds = sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
            if ds is None:
                raise ValueError(f"Dataset '{dataset_name}' not found (no current version).")

            inv = sess.query(self._InventoryItem).filter_by(name=inventory_name).first()
            if inv is None:
                raise ValueError(f"Inventory item '{inventory_name}' not found.")

            # Check for duplicate modality_key
            existing = (
                sess.query(self._DatasetMember)
                .filter_by(dataset_id=ds.id, modality_key=modality_key)
                .first()
            )
            if existing:
                raise ValueError(
                    f"Modality key '{modality_key}' already exists in dataset "
                    f"'{dataset_name}' v{ds.version}."
                )

            member = self._DatasetMember(
                dataset_id=ds.id,
                inventory_id=inv.id,
                modality_key=modality_key,
            )
            sess.add(member)
            sess.commit()
            sess.refresh(member)
            logger.info(
                "Added member '%s' (%s) to dataset '%s' v%d",
                inventory_name,
                modality_key,
                dataset_name,
                ds.version,
            )
            return member.id

    def remove_dataset_member(
        self,
        dataset_name: str,
        modality_key: str,
    ) -> None:
        """Remove a member from the current version of a dataset."""
        with self._session() as sess:
            ds = sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
            if ds is None:
                raise ValueError(f"Dataset '{dataset_name}' not found (no current version).")

            member = (
                sess.query(self._DatasetMember)
                .filter_by(dataset_id=ds.id, modality_key=modality_key)
                .first()
            )
            if member is None:
                raise ValueError(
                    f"Modality key '{modality_key}' not found in dataset "
                    f"'{dataset_name}' v{ds.version}."
                )
            sess.delete(member)
            sess.commit()
            logger.info(
                "Removed member '%s' from dataset '%s' v%d",
                modality_key,
                dataset_name,
                ds.version,
            )

    def get_dataset_members(
        self,
        dataset_name: str,
        version: int | None = None,
    ) -> list[dict[str, Any]]:
        """Return list of member dicts (including inventory item details)."""
        with self._session() as sess:
            if version is not None:
                ds = sess.query(self._Dataset).filter_by(name=dataset_name, version=version).first()
            else:
                ds = sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
                if ds is None:
                    ds = (
                        sess.query(self._Dataset)
                        .filter_by(name=dataset_name)
                        .order_by(self._Dataset.version.desc())
                        .first()
                    )
            if ds is None:
                return []

            members = sess.query(self._DatasetMember).filter_by(dataset_id=ds.id).all()
            result = []
            for m in members:
                inv = sess.get(self._InventoryItem, m.inventory_id)
                d = {
                    "member_id": m.id,
                    "dataset_id": m.dataset_id,
                    "modality_key": m.modality_key,
                    "inventory_id": m.inventory_id,
                }
                if inv:
                    d["inventory_name"] = inv.name
                    d["inventory_uri"] = inv.uri
                    d["inventory_modality"] = inv.modality
                    d["inventory_platform"] = inv.platform
                    d["inventory_format"] = inv.format
                    d["n_obs"] = inv.n_obs
                    d["n_vars"] = inv.n_vars
                result.append(d)
            return result

    def bump_dataset_version(self, dataset_name: str) -> int:
        """Create a new version of a dataset, copying all members. Returns new dataset id."""
        with self._session() as sess:
            current = (
                sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
            )
            if current is None:
                raise ValueError(f"Dataset '{dataset_name}' not found (no current version).")

            old_version = current.version
            current.is_current = False
            current.updated_at = _utcnow()

            now = _utcnow()
            new_ds = self._Dataset(
                name=dataset_name,
                version=old_version + 1,
                description=current.description,
                uri=current.uri,
                format=current.format,
                n_obs=current.n_obs,
                size_mb=current.size_mb,
                is_current=True,
                meta=dict(current.meta or {}),
                created_at=now,
                updated_at=now,
            )
            sess.add(new_ds)
            sess.flush()  # get new_ds.id

            # Copy members
            old_members = sess.query(self._DatasetMember).filter_by(dataset_id=current.id).all()
            for m in old_members:
                new_member = self._DatasetMember(
                    dataset_id=new_ds.id,
                    inventory_id=m.inventory_id,
                    modality_key=m.modality_key,
                )
                sess.add(new_member)

            sess.commit()
            sess.refresh(new_ds)
            logger.info(
                "Bumped dataset '%s' from v%d to v%d (id=%d)",
                dataset_name,
                old_version,
                new_ds.version,
                new_ds.id,
            )
            return new_ds.id

    # ------------------------------------------------------------------
    # Project-Dataset links (NEW)
    # ------------------------------------------------------------------

    def link_project_dataset(
        self,
        project_name: str,
        dataset_name: str,
        *,
        role: str = "primary",
        notes: str | None = None,
    ) -> int:
        """Link a project to a dataset. Returns link id. Idempotent."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found.")

            ds = sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
            if ds is None:
                raise ValueError(f"Dataset '{dataset_name}' not found (no current version).")

            existing = (
                sess.query(self._ProjectDataset)
                .filter_by(project_id=proj.id, dataset_id=ds.id)
                .first()
            )
            if existing:
                return existing.id

            link = self._ProjectDataset(
                project_id=proj.id,
                dataset_id=ds.id,
                role=role,
                notes=notes,
            )
            sess.add(link)
            sess.commit()
            sess.refresh(link)
            logger.info(
                "Linked project '%s' -> dataset '%s' (role=%s)",
                project_name,
                dataset_name,
                role,
            )
            return link.id

    def list_project_datasets(self, project_name: str) -> list[dict[str, Any]]:
        """Return all datasets linked to a project."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return []
            links = sess.query(self._ProjectDataset).filter_by(project_id=proj.id).all()
            result = []
            for lnk in links:
                ds = sess.get(self._Dataset, lnk.dataset_id)
                if ds:
                    d = self._to_dict(ds)
                    d["link_role"] = lnk.role
                    d["link_notes"] = lnk.notes
                    d["link_id"] = lnk.id
                    result.append(d)
            return result

    # ------------------------------------------------------------------
    # Project Phases (REWRITTEN -- uses project_phases table)
    # ------------------------------------------------------------------

    def upsert_phase(
        self,
        project_name: str,
        dataset_name: str,
        phase_group: str,
        subphase: str,
        *,
        status: str = "in_progress",
        uri: str | None = None,
        n_obs: int | None = None,
        n_vars: int | None = None,
        notes: str | None = None,
        metadata: dict | None = None,
    ) -> None:
        """Create or update a phase in project_phases.

        Parameters
        ----------
        project_name:
            Project name (must exist).
        dataset_name:
            Dataset name (must exist, uses is_current version).
        phase_group:
            Phase group: ``'data_processing'`` or ``'discovery'``.
        subphase:
            Subphase slug (e.g. ``'qc_filter'``, ``'clustering_v1'``).
        status:
            Phase status string.
        uri:
            Optional checkpoint URI.
        n_obs:
            Number of observations.
        n_vars:
            Number of variables.
        notes:
            Free-text notes (stored in metadata).
        metadata:
            Additional metadata dict.
        """
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")

            ds = sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
            if ds is None:
                raise ValueError(
                    f"Dataset '{dataset_name}' not found (no current version). "
                    "Create the dataset and link it to the project first."
                )

            existing = (
                sess.query(self._ProjectPhase)
                .filter_by(
                    project_id=proj.id,
                    dataset_id=ds.id,
                    phase_group=phase_group,
                    subphase=subphase,
                )
                .first()
            )

            now = _utcnow()
            if existing is not None:
                existing.status = status
                existing.updated_at = now
                if uri is not None:
                    existing.uri = uri
                if n_obs is not None:
                    existing.n_obs = n_obs
                if n_vars is not None:
                    existing.n_vars = n_vars
                meta = dict(existing.meta or {})
                if notes is not None:
                    meta["notes"] = notes
                if metadata:
                    meta.update(metadata)
                existing.meta = meta
            else:
                meta = dict(metadata or {})
                if notes is not None:
                    meta["notes"] = notes
                row = self._ProjectPhase(
                    project_id=proj.id,
                    dataset_id=ds.id,
                    phase_group=phase_group,
                    subphase=subphase,
                    status=status,
                    uri=uri,
                    n_obs=n_obs,
                    n_vars=n_vars,
                    meta=meta,
                    created_at=now,
                    updated_at=now,
                )
                sess.add(row)
            sess.commit()

    def get_phase(
        self,
        project_name: str,
        phase_group: str,
        subphase: str,
    ) -> dict[str, Any] | None:
        """Return phase info from project_phases."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return None
            row = (
                sess.query(self._ProjectPhase)
                .filter_by(
                    project_id=proj.id,
                    phase_group=phase_group,
                    subphase=subphase,
                )
                .first()
            )
            if row is None:
                return None
            meta = row.meta or {}
            return {
                "phase_group": row.phase_group,
                "subphase": row.subphase,
                "status": row.status,
                "uri": row.uri,
                "n_obs": row.n_obs,
                "n_vars": row.n_vars,
                "dataset_id": row.dataset_id,
                "notes": meta.get("notes"),
                "updated_at": row.updated_at,
            }

    def list_phases(
        self,
        project_name: str,
        phase_group: str | None = None,
    ) -> list[dict[str, Any]]:
        """Return phase summaries for a project from project_phases."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return []

            q = sess.query(self._ProjectPhase).filter_by(project_id=proj.id)
            if phase_group is not None:
                q = q.filter_by(phase_group=phase_group)

            rows = q.order_by(
                self._ProjectPhase.phase_group,
                self._ProjectPhase.subphase,
            ).all()

            result = []
            for r in rows:
                rmeta = r.meta or {}
                result.append(
                    {
                        "phase_group": r.phase_group,
                        "subphase": r.subphase,
                        "status": r.status,
                        "uri": r.uri,
                        "n_obs": r.n_obs,
                        "n_vars": r.n_vars,
                        "dataset_id": r.dataset_id,
                        "notes": rmeta.get("notes"),
                        "updated_at": r.updated_at,
                    }
                )
            return result

    def mark_phase_complete(
        self,
        project_name: str,
        phase_group: str,
        subphase: str,
    ) -> None:
        """Mark a phase as 'ready' in project_phases."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                raise ValueError(f"Project '{project_name}' not found in registry")

            rows = (
                sess.query(self._ProjectPhase)
                .filter_by(
                    project_id=proj.id,
                    phase_group=phase_group,
                    subphase=subphase,
                )
                .all()
            )
            if not rows:
                raise ValueError(
                    f"Phase '{phase_group}/{subphase}' not found for project '{project_name}'"
                )
            now = _utcnow()
            for r in rows:
                r.status = "ready"
                r.updated_at = now
            sess.commit()

    # ------------------------------------------------------------------
    # Provenance (NEW)
    # ------------------------------------------------------------------

    def record_provenance(
        self,
        tool: str,
        *,
        tool_version: str | None = None,
        reference_genome: str | None = None,
        reference_dataset: str | None = None,
        signature_source: str | None = None,
        n_input_obs: int | None = None,
        n_output_obs: int | None = None,
        params: dict | None = None,
        environment: dict | None = None,
        script_uri: str | None = None,
        agent: str | None = None,
        phase_id: int | None = None,
        dataset_id: int | None = None,
    ) -> int:
        """Record provenance for a transformation. Returns provenance id.

        Exactly one of phase_id or dataset_id must be provided.
        """
        # Validate exactly one FK
        fk_count = sum(x is not None for x in (phase_id, dataset_id))
        if fk_count != 1:
            raise ValueError(
                "Exactly one of phase_id or dataset_id must be provided. "
                f"Got {fk_count}."
            )

        if phase_id is not None:
            target_type = "phase"
        else:
            target_type = "dataset"

        with self._session() as sess:
            row = self._Provenance(
                target_type=target_type,
                phase_id=phase_id,
                dataset_id=dataset_id,
                tool=tool,
                tool_version=tool_version,
                reference_genome=reference_genome,
                reference_dataset=reference_dataset,
                signature_source=signature_source,
                n_input_obs=n_input_obs,
                n_output_obs=n_output_obs,
                params=params or {},
                environment=environment or {},
                script_uri=script_uri,
                agent=agent,
                created_at=_utcnow(),
            )
            sess.add(row)
            sess.commit()
            sess.refresh(row)
            logger.info(
                "Recorded provenance id=%d (tool=%s, target_type=%s)",
                row.id,
                tool,
                target_type,
            )
            return row.id

    def get_provenance(
        self,
        *,
        phase_id: int | None = None,
        dataset_id: int | None = None,
    ) -> list[dict[str, Any]]:
        """Return provenance records for a given target."""
        with self._session() as sess:
            q = sess.query(self._Provenance)
            if phase_id is not None:
                q = q.filter_by(target_type="phase", phase_id=phase_id)
            elif dataset_id is not None:
                q = q.filter_by(target_type="dataset", dataset_id=dataset_id)
            else:
                # No filter -- return all
                pass
            return [self._to_dict(r) for r in q.order_by(self._Provenance.id).all()]

    # ------------------------------------------------------------------
    # Data sources (uses new data_sources table)
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
        metadata: dict | None = None,
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
                meta=metadata or {},
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

    # ------------------------------------------------------------------
    # Datasets query (REWRITTEN)
    # ------------------------------------------------------------------

    def list_datasets(
        self,
        project_name: str | None = None,
    ) -> list[dict[str, Any]]:
        """List datasets, optionally filtered by project via project_datasets."""
        with self._session() as sess:
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                links = sess.query(self._ProjectDataset).filter_by(project_id=proj.id).all()
                dataset_ids = [lnk.dataset_id for lnk in links]
                if not dataset_ids:
                    return []
                rows = (
                    sess.query(self._Dataset)
                    .filter(self._Dataset.id.in_(dataset_ids))
                    .order_by(self._Dataset.name, self._Dataset.version)
                    .all()
                )
            else:
                rows = (
                    sess.query(self._Dataset)
                    .order_by(self._Dataset.name, self._Dataset.version)
                    .all()
                )
            return [self._to_dict(r) for r in rows]

    def get_dataset_uri(
        self,
        project_name: str,
        dataset_name: str,
    ) -> str | None:
        """Return the URI for a dataset linked to a project."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return None
            # Find the current version of the dataset
            ds = sess.query(self._Dataset).filter_by(name=dataset_name, is_current=True).first()
            if ds is None:
                return None
            # Verify it is linked to the project
            link = (
                sess.query(self._ProjectDataset)
                .filter_by(project_id=proj.id, dataset_id=ds.id)
                .first()
            )
            if link is None:
                return None
            return ds.uri

    def project_data_summary(self, project_name: str) -> dict[str, Any]:
        """Return summary of datasets and their members for a project."""
        with self._session() as sess:
            proj = sess.query(self._Project).filter_by(name=project_name).first()
            if proj is None:
                return {
                    "project": project_name,
                    "n_datasets": 0,
                    "n_inventory_items": 0,
                    "datasets": [],
                }

            links = sess.query(self._ProjectDataset).filter_by(project_id=proj.id).all()
            datasets_info = []
            all_inventory_ids: set[int] = set()

            for lnk in links:
                ds = sess.get(self._Dataset, lnk.dataset_id)
                if ds is None:
                    continue
                members = sess.query(self._DatasetMember).filter_by(dataset_id=ds.id).all()
                member_info = []
                for m in members:
                    all_inventory_ids.add(m.inventory_id)
                    inv = sess.get(self._InventoryItem, m.inventory_id)
                    member_info.append(
                        {
                            "modality_key": m.modality_key,
                            "inventory_name": inv.name if inv else None,
                            "modality": inv.modality if inv else None,
                            "n_obs": inv.n_obs if inv else None,
                        }
                    )
                datasets_info.append(
                    {
                        "name": ds.name,
                        "version": ds.version,
                        "is_current": ds.is_current,
                        "format": ds.format,
                        "role": lnk.role,
                        "n_members": len(members),
                        "members": member_info,
                    }
                )

            return {
                "project": project_name,
                "n_datasets": len(datasets_info),
                "n_inventory_items": len(all_inventory_ids),
                "datasets": datasets_info,
            }

    # ------------------------------------------------------------------
    # Patients (UNCHANGED)
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
                q = q.filter(self._Patient.meta["diagnosis"].astext.ilike(f"%{diagnosis}%"))
            if tissue is not None:
                q = q.filter(self._Patient.meta["tissue_of_origin"].astext.ilike(f"%{tissue}%"))
            return [self._to_dict(r) for r in q.order_by(self._Patient.patient_id).all()]

    def list_subjects(
        self,
        project_name: str | None = None,
        diagnosis: str | None = None,
        tissue: str | None = None,
    ) -> list[dict[str, Any]]:
        """Backward-compatible alias for :meth:`list_patients`.

        The project_name filter checks for patients linked via
        patient -> sample -> inventory_item -> dataset_member -> dataset -> project.
        """
        with self._session() as sess:
            q = sess.query(self._Patient)
            if project_name is not None:
                proj = sess.query(self._Project).filter_by(name=project_name).first()
                if proj is None:
                    return []
                # patient -> sample -> inventory_item -> dataset_member -> project_dataset
                patient_ids_subq = (
                    sess.query(self._Sample.patient_id)
                    .join(self._InventoryItem, self._InventoryItem.sample_id == self._Sample.id)
                    .join(self._DatasetMember, self._DatasetMember.inventory_id == self._InventoryItem.id)
                    .join(self._Dataset, self._DatasetMember.dataset_id == self._Dataset.id)
                    .join(self._ProjectDataset, self._ProjectDataset.dataset_id == self._Dataset.id)
                    .filter(self._ProjectDataset.project_id == proj.id)
                    .filter(self._Sample.patient_id.isnot(None))
                    .distinct()
                    .scalar_subquery()
                )
                q = q.filter(self._Patient.id.in_(patient_ids_subq))
            if diagnosis is not None:
                q = q.filter(self._Patient.meta["diagnosis"].astext.ilike(f"%{diagnosis}%"))
            if tissue is not None:
                q = q.filter(self._Patient.meta["tissue_of_origin"].astext.ilike(f"%{tissue}%"))
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
    # Samples
    # ------------------------------------------------------------------

    def register_sample(
        self,
        sample_id: str,
        *,
        patient_id: str | None = None,
        tissue: str | None = None,
        collection_date: str | None = None,
        metadata: dict | None = None,
    ) -> int:
        """Register a sample. Returns sample DB id. Idempotent by sample_id."""
        with self._session() as sess:
            existing = sess.query(self._Sample).filter_by(sample_id=sample_id).first()
            if existing:
                logger.info("Sample '%s' already registered (id=%d)", sample_id, existing.id)
                return existing.id

            pat_db_id = None
            if patient_id is not None:
                pat = sess.query(self._Patient).filter_by(patient_id=patient_id).first()
                if pat is None:
                    raise ValueError(f"Patient '{patient_id}' not found.")
                pat_db_id = pat.id

            sample = self._Sample(
                sample_id=sample_id,
                patient_id=pat_db_id,
                tissue=tissue,
                collection_date=collection_date,
                meta=metadata or {},
            )
            sess.add(sample)
            sess.commit()
            sess.refresh(sample)
            logger.info("Registered sample '%s' (id=%d)", sample_id, sample.id)
            return sample.id

    def get_sample(self, sample_id: str) -> dict[str, Any] | None:
        """Return sample dict by sample_id, or None."""
        with self._session() as sess:
            row = sess.query(self._Sample).filter_by(sample_id=sample_id).first()
            return self._to_dict(row) if row else None

    def list_samples(
        self,
        patient_id: str | None = None,
        tissue: str | None = None,
    ) -> list[dict[str, Any]]:
        """List samples with optional filters."""
        with self._session() as sess:
            q = sess.query(self._Sample)
            if patient_id is not None:
                pat = sess.query(self._Patient).filter_by(patient_id=patient_id).first()
                if pat is None:
                    return []
                q = q.filter_by(patient_id=pat.id)
            if tissue is not None:
                q = q.filter(self._Sample.tissue.ilike(f"%{tissue}%"))
            return [self._to_dict(r) for r in q.order_by(self._Sample.sample_id).all()]

    # ------------------------------------------------------------------
    # Status summary (REWRITTEN)
    # ------------------------------------------------------------------

    def status(self) -> dict[str, Any]:
        """Return a high-level summary of current registry state."""
        with self._session() as sess:
            n_projects = sess.query(self._Project).count()
            n_data_sources = sess.query(self._DataSource).count()
            n_inventory_items = sess.query(self._InventoryItem).count()
            n_datasets = sess.query(self._Dataset).filter_by(is_current=True).count()
            n_patients = sess.query(self._Patient).count()
            n_samples = sess.query(self._Sample).count()
            projects = [p.name for p in sess.query(self._Project).filter_by(status="active").all()]

            # Per-project phase summary from project_phases
            from sqlalchemy import func

            phase_summary: dict[str, dict[str, int]] = {}
            for proj_name in projects:
                proj = sess.query(self._Project).filter_by(name=proj_name).first()
                if proj is None:
                    continue
                phase_counts = (
                    sess.query(self._ProjectPhase.status, func.count())
                    .filter_by(project_id=proj.id)
                    .group_by(self._ProjectPhase.status)
                    .all()
                )
                counts = dict(phase_counts)
                if counts:
                    phase_summary[proj_name] = counts

        return {
            "n_projects": n_projects,
            "n_data_sources": n_data_sources,
            "n_inventory_items": n_inventory_items,
            "n_datasets": n_datasets,
            "n_patients": n_patients,
            "n_samples": n_samples,
            # Backward-compat keys
            "n_subjects": n_patients,
            "n_data": n_inventory_items,
            "n_biodata": n_inventory_items,
            "active_slurm_jobs": 0,
            "running_agent_tasks": 0,
            "active_projects": projects,
            "phase_summary": phase_summary,
        }

    # ------------------------------------------------------------------
    # Backward compatibility wrappers (DEPRECATED)
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
        """Deprecated. Use register_inventory_item() instead."""
        warnings.warn(
            "register_data() is deprecated. Use register_inventory_item() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        name = f"{project_name}_{phase}_{_utcnow().replace(':', '-')}"
        return self.register_inventory_item(
            name=name,
            uri=uri,
            modality=category or "unknown",
            platform=platform,
            fmt=fmt or "h5ad",
            n_obs=n_obs,
            n_vars=n_vars,
            size_mb=size_mb,
            metadata=metadata,
        )

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
        """Deprecated. Use register_inventory_item() instead."""
        warnings.warn(
            "register_dataset() is deprecated. Use register_inventory_item() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        meta: dict[str, Any] = {}
        if md5:
            meta["md5"] = md5
        if validated:
            meta["validated"] = True
        name = f"{project_name}_{phase}_{_utcnow().replace(':', '-')}"
        return self.register_inventory_item(
            name=name,
            uri=uri,
            modality="unknown",
            fmt=fmt,
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
        """Deprecated. Use register_inventory_item() instead."""
        warnings.warn(
            "register_biodata() is deprecated. Use register_inventory_item() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        meta: dict[str, Any] = {}
        if md5:
            meta["md5"] = md5
        if type_kwargs:
            meta.update(type_kwargs)
        name = f"{project_name}_{category}_{platform}_{_utcnow().replace(':', '-')}"
        return self.register_inventory_item(
            name=name,
            uri=uri,
            modality=category,
            platform=platform,
            fmt=fmt or "h5ad",
            n_obs=n_obs,
            n_vars=n_vars,
            size_mb=size_mb,
            metadata=meta if meta else None,
        )

    def list_biodata(
        self,
        project_name: str | None = None,
        category: str | None = None,
        platform: str | None = None,
        sample_db_id: int | None = None,
        phase: str | None = None,
    ) -> list[dict[str, Any]]:
        """Deprecated. Use list_inventory_items() instead."""
        warnings.warn(
            "list_biodata() is deprecated. Use list_inventory_items() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.list_inventory_items(
            modality=category,
            platform=platform,
        )

    def link_project_data_source(
        self,
        project_name: str,
        data_source_name: str,
        *,
        role: str = "input",
        notes: str | None = None,
    ) -> int:
        """Deprecated. Use link_project_dataset() instead."""
        warnings.warn(
            "link_project_data_source() is deprecated. Use link_project_dataset() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.link_project_dataset(
            project_name=project_name,
            dataset_name=data_source_name,
            role=role,
            notes=notes,
        )

    def list_project_data_sources(self, project_name: str) -> list[dict[str, Any]]:
        """Deprecated. Use list_project_datasets() instead."""
        warnings.warn(
            "list_project_data_sources() is deprecated. Use list_project_datasets() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.list_project_datasets(project_name)

    # ------------------------------------------------------------------
    # Legacy methods kept for burn-in period
    # ------------------------------------------------------------------

    def get_biodata(self, biodata_id: int) -> dict[str, Any] | None:
        """Return an inventory item dict by id."""
        with self._session() as sess:
            row = sess.get(self._InventoryItem, biodata_id)
            if row is None:
                return None
            return self._to_dict(row)

    def update_dataset_status(self, dataset_id: int, status: str) -> None:
        """Update dataset status."""
        with self._session() as sess:
            row = sess.get(self._Dataset, dataset_id)
            if row is None:
                raise ValueError(f"Dataset id={dataset_id} not found")
            row.updated_at = _utcnow()
            sess.commit()

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
    print(f"  Data sources      : {s['n_data_sources']}")
    print(f"  Inventory items   : {s['n_inventory_items']}")
    print(f"  Datasets          : {s['n_datasets']}")
    print(f"  Patients          : {s['n_patients']}")
    print(f"  Samples           : {s['n_samples']}")
    if s["active_projects"]:
        print("\n  Active projects:")
        for proj_name in s["active_projects"]:
            print(f"    - {proj_name}")
            phase_summary = s.get("phase_summary", {}).get(proj_name, {})
            if phase_summary:
                parts = ", ".join(f"{k}={v}" for k, v in sorted(phase_summary.items()))
                print(f"      phase status: {parts}")
            datasets = reg.list_project_datasets(project_name=proj_name)
            for ds in datasets[-5:]:
                role = ds.get("link_role", "primary")
                print(
                    f"        dataset={ds['name']} v{ds['version']} role={role} "
                    f"format={ds['format']}"
                )
    print()


__all__ = ["Registry"]
