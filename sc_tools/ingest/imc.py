"""IMC pipeline command builder and image utilities for Phase 0a/0b.

Builds shell commands to run the ElementoLab IMC pipeline for segmentation
and single-cell extraction from MCD files. The pipeline outputs to::

    processed/{sample}/tiffs/   # multi-channel TIFF stacks and channel CSVs
    processed/{sample}/masks/   # segmentation masks
    processed/{sample}/cells.h5ad  # single-cell quantification (Phase 0b input)

**TIFF stack naming convention** (ElementoLab/steinbock)::

    {roi_id}_full.tiff          # (C, H, W) multi-channel intensity stack
    {roi_id}_full.csv           # channel index → MarkerName(IsotopeTag)
    {roi_id}_full_mask.tiff     # labeled cell segmentation mask
    {roi_id}_full_nucmask.tiff  # labeled nuclear segmentation mask
    {roi_id}_Probabilities.tiff # ilastik probability maps

**Channel CSV format** (``*_full.csv``)::

    ,channel
    0,CD3(Er170)
    1,PanCK(Pt195)
    ...

The channel name format is ``MarkerName(IsotopeTag)`` where ``MarkerName``
(before ``(``) is the protein target and ``IsotopeTag`` (inside ``()``) is
the metal isotope label.

**Panel/metadata CSV format** (``channel_labels.csv``)::

    channel,Target,Metal_Tag,Atom,full,ilastik
    CD3(Er170),CD3,Er170,170,1,1
    PanCK(Pt195),PanCK,Pt195,195,1,0
    ...

Phase 0b loading is handled by ``sc_tools.ingest.loaders.load_imc_sample()``,
which reads ``processed/{sample}/cells.h5ad`` and optionally the TIFF stacks.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def build_imc_pipeline_cmd(
    mcd_file: str,
    panel_csv: str,
    output_dir: str,
    *,
    pipeline_dir: str | None = None,
) -> str:
    """Build shell command to run the ElementoLab IMC pipeline.

    Parameters
    ----------
    mcd_file
        Path to .mcd file.
    panel_csv
        Path to panel CSV defining channels.
    output_dir
        Root output directory. The pipeline writes per-sample results under
        ``output_dir/{sample}/tiffs/``, ``masks/``, and ``cells.h5ad``.
        Typically set to ``processed/`` in the project directory.
    pipeline_dir
        Path to cloned ElementoLab IMC pipeline repo. Defaults to
        ``~/elementolab/imc``.

    Returns
    -------
    Shell command string.
    """
    if pipeline_dir is None:
        pipeline_dir = "~/elementolab/imc"

    parts = [
        f"python {pipeline_dir}/run_pipeline.py",
        f"--mcd {mcd_file}",
        f"--panel {panel_csv}",
        f"--output {output_dir}",
    ]
    return " \\\n    ".join(parts)


def _parse_channel_name(channel_str: str) -> tuple[str, str]:
    """Parse ``MarkerName(IsotopeTag)`` into ``(protein_name, isotope_tag)``.

    Examples
    --------
    >>> _parse_channel_name("CD3(Er170)")
    ('CD3', 'Er170')
    >>> _parse_channel_name("DNA1(Ir191)")
    ('DNA1', 'Ir191')
    >>> _parse_channel_name("<EMPTY>(In115)")
    ('<EMPTY>', 'In115')
    >>> _parse_channel_name("CD3")
    ('CD3', '')
    """
    m = re.match(r"^(.+?)\((.+)\)$", channel_str.strip())
    if m:
        return m.group(1).strip(), m.group(2).strip()
    return channel_str.strip(), ""


class IMCPanelMapper:
    """Map protein/marker names to channel indices in an IMC TIFF stack.

    IMC pipelines (ElementoLab/steinbock) produce a per-ROI multi-channel TIFF
    stack (``*_full.tiff``) and a companion CSV (``*_full.csv``) that maps each
    channel index to a string like ``CD3(Er170)`` or ``PanCK(Pt195)``.

    This class parses the channel CSV and resolves arbitrary name queries
    (protein name, isotope tag, full string, or partial name) to a channel
    index for slicing the TIFF stack.

    Channel name parsing: ``MarkerName(IsotopeTag)``
    - protein name = text before ``(``
    - isotope tag  = text inside ``()``

    Matching precedence:
    1. Exact match on protein name (Target column)
    2. Exact match on full channel string (``CD3(Er170)``)
    3. Case-insensitive exact match on protein name
    4. Exact match on isotope tag
    5. Partial / substring match on protein name
    6. ``None`` with a ``logger.warning()`` — never raises

    Parameters
    ----------
    full_csv
        Optional path to per-ROI channel index CSV (``*_full.csv``). If
        provided, used as the authoritative channel list for the TIFF stack.
    panel_csv
        Optional path to panel metadata CSV (``channel_labels.csv``) with
        columns ``channel, Target, Metal_Tag, Atom, full, ilastik``. Provides
        additional aliases and the ``ilastik`` / ``full`` flags.

    Examples
    --------
    >>> mapper = IMCPanelMapper(full_csv="ROI01_full.csv")
    >>> mapper.resolve("CD3")          # returns 33 (index in TIFF stack)
    33
    >>> mapper.build_rgb_indices()     # returns {R: 46, G: 33, B: 42}
    {'R': 46, 'G': 33, 'B': 42}
    """

    def __init__(
        self,
        full_csv: str | Path | None = None,
        panel_csv: str | Path | None = None,
    ) -> None:
        # Ordered list of channel records (preserves TIFF stack order)
        # Each record: {index, full_string, protein_name, isotope_tag}
        self._channels: list[dict] = []

        # Lookup maps: lower-case query string -> channel index
        self._protein_to_idx: dict[str, int] = {}
        self._full_str_to_idx: dict[str, int] = {}
        self._isotope_to_idx: dict[str, int] = {}

        # Additional flags from panel CSV (protein_name -> flag)
        self._ilastik_flags: dict[str, int] = {}
        self._full_flags: dict[str, int] = {}

        if full_csv is not None:
            self.from_full_csv(full_csv)
        if panel_csv is not None:
            self.from_panel_csv(panel_csv)

    # ------------------------------------------------------------------
    # Population helpers
    # ------------------------------------------------------------------

    def _register_channel(self, idx: int, full_string: str) -> None:
        """Register one channel into the lookup tables."""
        protein, isotope = _parse_channel_name(full_string)
        record = {
            "index": idx,
            "full_string": full_string,
            "protein_name": protein,
            "isotope_tag": isotope,
        }
        # Extend list if needed
        while len(self._channels) <= idx:
            self._channels.append({})
        self._channels[idx] = record

        self._full_str_to_idx[full_string.lower()] = idx
        if protein and protein != "<EMPTY>":
            self._protein_to_idx[protein.lower()] = idx
        if isotope:
            self._isotope_to_idx[isotope.lower()] = idx

    def from_full_csv(self, path: str | Path) -> IMCPanelMapper:
        """Parse a per-ROI ``*_full.csv`` channel index file.

        Expected format (two columns)::

            ,channel
            0,CD3(Er170)
            1,PanCK(Pt195)
            ...

        Parameters
        ----------
        path
            Path to the ``*_full.csv`` file.

        Returns
        -------
        self (for chaining)
        """
        import pandas as pd

        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Channel CSV not found: {path}")

        df = pd.read_csv(path, index_col=0)
        # Column could be named 'channel' or unnamed
        channel_col = df.columns[0]

        for idx_val, row in df.iterrows():
            try:
                idx = int(idx_val)
            except (ValueError, TypeError):
                continue
            self._register_channel(idx, str(row[channel_col]).strip())

        return self

    def from_panel_csv(self, path: str | Path) -> IMCPanelMapper:
        """Parse a panel metadata CSV (``channel_labels.csv``).

        Expected columns: ``channel, Target, Metal_Tag, Atom, full, ilastik``.
        The ``channel`` column contains the full string like ``CD3(Er170)``.
        ``Target`` is the protein name. Adds isotope-tag aliases and stores
        ``full`` / ``ilastik`` flags.

        If ``from_full_csv`` has not been called, this also registers channel
        order (0-based row index).

        Parameters
        ----------
        path
            Path to the panel CSV file.

        Returns
        -------
        self (for chaining)
        """
        import pandas as pd

        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Panel CSV not found: {path}")

        df = pd.read_csv(path)

        # Drop mcd_file column if present (aggregated CSV)
        if "mcd_file" in df.columns:
            df = df.drop(columns=["mcd_file"])

        channel_col = next((c for c in ["channel", "Channel"] if c in df.columns), None)
        target_col = next((c for c in ["Target", "target", "marker"] if c in df.columns), None)
        ilastik_col = next((c for c in ["ilastik"] if c in df.columns), None)
        full_col_name = next((c for c in ["full"] if c in df.columns), None)

        for row_idx, row in df.iterrows():
            ch_str = str(row[channel_col]).strip() if channel_col else ""
            protein = str(row[target_col]).strip() if target_col else ""

            if not protein or protein.lower() in ("nan", "none", "<empty>"):
                continue

            # Store flags
            if ilastik_col:
                self._ilastik_flags[protein.lower()] = int(row[ilastik_col])
            if full_col_name:
                self._full_flags[protein.lower()] = int(row[full_col_name])

            # If this channel is already registered via from_full_csv, add
            # protein-name alias. Otherwise register at row_idx.
            ch_lower = ch_str.lower()
            if ch_lower in self._full_str_to_idx:
                idx = self._full_str_to_idx[ch_lower]
                if protein.lower() not in self._protein_to_idx:
                    self._protein_to_idx[protein.lower()] = idx
            else:
                # Register by row order (panel CSV order = TIFF stack order)
                self._register_channel(int(row_idx), ch_str)
                if protein.lower() not in self._protein_to_idx:
                    _, isotope = _parse_channel_name(ch_str)
                    self._protein_to_idx[protein.lower()] = int(row_idx)

        return self

    def set_from_var_names(self, var_names: list[str]) -> IMCPanelMapper:
        """Register channel names from ``adata.var_names`` (no index CSV needed).

        Each var_name may be a full string like ``CD3(Er170)`` or a plain
        protein name. Assigns indices in list order (0, 1, 2, ...).

        Parameters
        ----------
        var_names
            Ordered list of marker names from ``adata.var_names``.

        Returns
        -------
        self (for chaining)
        """
        for idx, name in enumerate(var_names):
            self._register_channel(idx, str(name).strip())
        return self

    # ------------------------------------------------------------------
    # Resolution
    # ------------------------------------------------------------------

    def resolve(self, name: str) -> int | None:
        """Map any name form to a TIFF stack channel index.

        Parameters
        ----------
        name
            Query: protein name, isotope tag, full channel string, or partial.

        Returns
        -------
        Integer channel index (0-based) or ``None`` if unresolvable (logs a
        warning). Never raises.
        """
        if not name or name.lower() in ("<empty>", "nan", "none"):
            return None

        lo = name.lower()

        # 1. Exact protein name
        if lo in self._protein_to_idx:
            return self._protein_to_idx[lo]

        # 2. Exact full string
        if lo in self._full_str_to_idx:
            return self._full_str_to_idx[lo]

        # 3. Exact isotope tag
        if lo in self._isotope_to_idx:
            return self._isotope_to_idx[lo]

        # 4. Partial / substring match on protein names
        matches = [(k, v) for k, v in self._protein_to_idx.items() if lo in k or k in lo]
        if len(matches) == 1:
            return matches[0][1]
        if len(matches) > 1:
            # Prefer prefix match
            prefix = [(k, v) for k, v in matches if k.startswith(lo)]
            if len(prefix) == 1:
                return prefix[0][1]
            logger.warning(
                "IMCPanelMapper: ambiguous match for %r -> %s; returning first",
                name,
                [k for k, _ in matches],
            )
            return matches[0][1]

        logger.warning("IMCPanelMapper: could not resolve channel name %r", name)
        return None

    def resolve_name(self, name: str) -> str | None:
        """Return the canonical protein name for a query (instead of index)."""
        idx = self.resolve(name)
        if idx is None:
            return None
        ch = self._channels[idx] if idx < len(self._channels) else {}
        return ch.get("protein_name") or ch.get("full_string")

    def build_rgb_indices(
        self,
        r: str = "PanCK",
        g: str = "CD3",
        b: str = "DNA1",
    ) -> dict[str, int | None]:
        """Return TIFF stack indices for R/G/B composite channels.

        Parameters
        ----------
        r, g, b
            Requested marker names for red, green, blue channels.

        Returns
        -------
        dict with keys ``'R'``, ``'G'``, ``'B'`` and integer stack indices
        (or ``None`` if not found — logs a warning).
        """
        return {
            "R": self.resolve(r),
            "G": self.resolve(g),
            "B": self.resolve(b),
        }

    def channel_names(self) -> list[str]:
        """Return ordered protein names for all registered channels."""
        return [ch.get("protein_name") or ch.get("full_string", "") for ch in self._channels if ch]

    def channel_strings(self) -> list[str]:
        """Return ordered full channel strings (``MarkerName(IsotopeTag)``)."""
        return [ch.get("full_string", "") for ch in self._channels if ch]

    def n_channels(self) -> int:
        """Number of registered channels."""
        return len(self._channels)

    def list_available(self) -> list[str]:
        """Return all available protein names (excluding empty channels)."""
        return [
            ch["protein_name"]
            for ch in self._channels
            if ch and ch.get("protein_name") and ch["protein_name"] != "<EMPTY>"
        ]

    def get_ilastik_indices(self) -> list[int]:
        """Return indices of channels with ``ilastik=1`` (from panel CSV)."""
        return [
            self._protein_to_idx[k]
            for k in self._ilastik_flags
            if self._ilastik_flags[k] == 1 and k in self._protein_to_idx
        ]


# ---------------------------------------------------------------------------
# Image composite builder
# ---------------------------------------------------------------------------


def build_imc_composite(
    tiff_path: str | Path,
    channel_csv: str | Path,
    panel_mapper: IMCPanelMapper | None = None,
    *,
    r: str = "PanCK",
    g: str = "CD3",
    b: str = "DNA1",
    percentile: float = 99.5,
    downsample: int = 1,
    load_mask: bool = False,
    load_probabilities: bool = False,
) -> dict:
    """Load an IMC ``*_full.tiff`` stack and build a standard spatial image dict.

    Reads the multi-channel TIFF stack (``*_full.tiff``) using the companion
    channel CSV (``*_full.csv``) for channel name resolution, builds an
    arcsinh-normalized ``(C, H, W)`` full stack and a ``(H, W, 3)`` uint8 RGB
    composite from three user-specified channels.

    The returned dict matches the ``adata.uns['spatial'][library_id]`` schema
    expected by squidpy and ``sc.pl.spatial``.

    Parameters
    ----------
    tiff_path
        Path to the ``*_full.tiff`` multi-channel stack file.
    channel_csv
        Path to the companion ``*_full.csv`` channel index file.
    panel_mapper
        Pre-built ``IMCPanelMapper``. If ``None``, one is built from
        ``channel_csv`` automatically.
    r, g, b
        Marker names for the red, green, blue composite channels.
    percentile
        Percentile used for clipping each RGB channel before uint8 scaling.
    downsample
        Spatial downsampling factor (1 = no downsampling). Applied uniformly
        to H and W dimensions.
    load_mask
        If ``True``, also load ``*_full_mask.tiff`` from the same directory
        and store under ``images['mask']``.
    load_probabilities
        If ``True``, also load ``*_Probabilities.tiff`` and store under
        ``images['probabilities']``.

    Returns
    -------
    dict with keys ``'images'``, ``'scalefactors'``, and ``'metadata'``.

    Notes
    -----
    arcsinh normalization uses ``arcsinh(x / 5)`` per channel, matching the
    ``sc_tools.pp.normalize.arcsinh_transform`` convention.
    """
    try:
        import tifffile
    except ImportError as e:
        raise ImportError("tifffile is required for IMC image loading: pip install tifffile") from e

    tiff_path = Path(tiff_path)
    channel_csv = Path(channel_csv)

    if not tiff_path.exists():
        raise FileNotFoundError(f"TIFF stack not found: {tiff_path}")
    if not channel_csv.exists():
        raise FileNotFoundError(f"Channel CSV not found: {channel_csv}")

    # Build mapper if not provided
    if panel_mapper is None:
        panel_mapper = IMCPanelMapper(full_csv=channel_csv)
    else:
        # Ensure the CSV is loaded even if mapper was pre-built
        if panel_mapper.n_channels() == 0:
            panel_mapper.from_full_csv(channel_csv)

    # Load full TIFF stack: expected (C, H, W) or (H, W) for single channel
    raw = tifffile.imread(str(tiff_path)).astype(np.float32)
    if raw.ndim == 2:
        raw = raw[np.newaxis, ...]  # (1, H, W)
    elif raw.ndim != 3:
        raise ValueError(f"Unexpected TIFF shape {raw.shape}; expected (C, H, W)")

    if downsample > 1:
        raw = raw[:, ::downsample, ::downsample]

    # arcsinh-normalize full stack
    full_norm = np.arcsinh(raw / 5.0)

    # Build RGB composite
    rgb_indices = panel_mapper.build_rgb_indices(r=r, g=g, b=b)
    ch_names = panel_mapper.channel_names()

    def _percentile_scale(idx: int | None) -> np.ndarray:
        h, w = full_norm.shape[1], full_norm.shape[2]
        if idx is None or idx >= full_norm.shape[0]:
            return np.zeros((h, w), dtype=np.uint8)
        arr = full_norm[idx]
        p_high = float(np.percentile(arr, percentile))
        if p_high == 0:
            return np.zeros((h, w), dtype=np.uint8)
        return np.clip((arr / p_high) * 255, 0, 255).astype(np.uint8)

    hires = np.stack(
        [
            _percentile_scale(rgb_indices.get("R")),
            _percentile_scale(rgb_indices.get("G")),
            _percentile_scale(rgb_indices.get("B")),
        ],
        axis=-1,
    )  # (H, W, 3) uint8

    scalef = 1.0 / downsample
    images = {"hires": hires, "full": full_norm}

    # Optionally load mask and probabilities from same directory
    tiff_dir = tiff_path.parent
    # Infer ROI stem from TIFF filename (remove _full.tiff suffix)
    stem = tiff_path.name
    if stem.endswith("_full.tiff"):
        roi_stem = stem[: -len("_full.tiff")]
    elif stem.endswith("_full.tif"):
        roi_stem = stem[: -len("_full.tif")]
    else:
        roi_stem = tiff_path.stem

    if load_mask:
        for mask_suffix in ["_full_mask.tiff", "_full_mask.tif"]:
            mask_path = tiff_dir / f"{roi_stem}{mask_suffix}"
            if mask_path.exists():
                try:
                    images["mask"] = tifffile.imread(str(mask_path))
                except Exception as exc:
                    logger.warning("Could not load mask %s: %s", mask_path, exc)
                break
        else:
            logger.warning("Cell mask not found for ROI %s", roi_stem)

    if load_probabilities:
        for prob_suffix in ["_Probabilities.tiff", "_Probabilities.tif"]:
            prob_path = tiff_dir / f"{roi_stem}{prob_suffix}"
            if prob_path.exists():
                try:
                    images["probabilities"] = tifffile.imread(str(prob_path)).astype(np.float32)
                except Exception as exc:
                    logger.warning("Could not load probabilities %s: %s", prob_path, exc)
                break
        else:
            logger.warning("Probability map not found for ROI %s", roi_stem)

    return {
        "images": images,
        "scalefactors": {
            "tissue_hires_scalef": scalef,
            "tissue_lowres_scalef": scalef,
            "spot_diameter_fullres": 1.0,  # IMC: 1 pixel ~ 1 µm
        },
        "metadata": {
            "channels": ch_names,
            "channel_strings": panel_mapper.channel_strings(),
            "rgb_channels": {
                "R": ch_names[rgb_indices["R"]]
                if rgb_indices.get("R") is not None and rgb_indices["R"] < len(ch_names)
                else None,
                "G": ch_names[rgb_indices["G"]]
                if rgb_indices.get("G") is not None and rgb_indices["G"] < len(ch_names)
                else None,
                "B": ch_names[rgb_indices["B"]]
                if rgb_indices.get("B") is not None and rgb_indices["B"] < len(ch_names)
                else None,
            },
            "rgb_indices": rgb_indices,
            "pixel_size_um": 1.0,
        },
    }
