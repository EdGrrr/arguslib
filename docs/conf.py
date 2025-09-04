# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
from pathlib import Path
import sys

# Make the project root directory available to Sphinx.
# This allows autodoc to find the `arguslib` package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "arguslib"
copyright = "2025, Oliver G. A. Driver, Edward Gryspeerdt"
author = "Oliver G. A. Driver, Edward Gryspeerdt"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",
    "sphinx_gallery.gen_gallery",
    "matplotlib.sphinxext.roles",
]

# Tell Sphinx to mock heavy or problematic imports so that the documentation
# can be built without having to install them all. This is crucial for
# services like Read the Docs and for avoiding dependency installation issues.
# NOTE: For sphinx-gallery to run examples, dependencies must be installed, not mocked.
autodoc_mock_imports = [
    # "matplotlib",
    # "numpy",
    # "pyart",
    # "yaml",
    # "xarray",
    # "cv2",
    # "netCDF4",
    # "tqdm",
    # "csat2",
    # "pandas",
    # "PIL",
    # "pytz",
    # "typing_extensions",
    # "sphinx.ext.autodoc.canonical_imports",
]

# -- Sphinx-gallery settings -------------------------------------------------
from sphinx_gallery.sorting import FileNameSortKey

sphinx_gallery_conf = {
    "examples_dirs": "examples",  # path to your example scripts
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated rst files
    "within_subsection_order": FileNameSortKey,
    "filename_pattern": r"\.py$",
    "ignore_pattern": r"__init__\.py$",
    "doc_module": ("arguslib",),
}

# -- Autodoc settings --------------------------------------------------------
autodoc_member_order = "bysource"
autosummary_generate = True
autodoc_default_options = {
    # 'members': True,
    "undoc-members": True,
    "show-inheritance": True,
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]

html_theme_options = {
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink"],
}
