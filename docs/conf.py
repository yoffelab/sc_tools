"""Sphinx configuration for sc_tools documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- Project info -------------------------------------------------------------
project = "sc_tools"
author = "Junbum Kim"
release = "0.1.0"
copyright = "2024-2026, Junbum Kim"  # noqa: A001

# -- Extensions ---------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_design",
    "myst_nb",
]

# -- Autodoc ------------------------------------------------------------------
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "private-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"

# autodoc_mock_imports — keep this list minimal.
#
# All of sc_tools's optional soft dependencies (scvi, tangram, gseapy, torch,
# rapids_singlecell, cupy …) use try/except ImportError guards, so they fail
# gracefully when not installed.  Mocking packages that are absent can cause
# secondary failures: e.g. mocking torch sets torch.__spec__ = None, which
# makes anndata's find_spec("torch") raise ValueError; mocking cupy triggers
# anndata's is_cupy_importable() → True → singledispatch TypeError.
#
# Only mock packages that ARE installed on the build machine but must not be
# imported during docs generation (e.g. packages requiring a GPU at import
# time).  On a standard CPU-only docs environment this list can be empty.
autodoc_mock_imports: list[str] = []

# -- Napoleon (NumPy-style only) ----------------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# -- myst-nb (static notebooks, no execution) ---------------------------------
nb_execution_mode = "off"
myst_enable_extensions = ["colon_fence"]

# -- Intersphinx --------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", None),
}

# -- Exclude patterns (suppress myst-nb jupyter_execute toctree warnings) ----
exclude_patterns = ["_build", "_build/jupyter_execute/**"]

# -- HTML theme ---------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "github_url": "https://github.com/yoffelab/sc_tools",
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "secondary_sidebar_items": ["page-toc", "sourcelink"],
    "show_nav_level": 2,
    "navigation_depth": 3,
    "icon_links": [],
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "sc_tools"
html_show_sourcelink = True
