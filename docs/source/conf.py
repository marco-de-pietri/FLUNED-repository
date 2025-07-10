# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import sys
import importlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]   
sys.path.insert(0, str(ROOT / "src"))

pkg_name = "fluned"                      # src/modulename/__init__.py
pkg = importlib.import_module(pkg_name)

release = version = pkg.__version__          # e.g., "1.4.3"

# sys.path.insert(
#     0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
# )
#
# try:
#     from importlib.metadata import version as _get_version
# except ImportError:
#     from pkg_resources import get_distribution as _get_version
#
# release = _get_version("fluned")
# # short X.Y version
# version = ".".join(release.split(".")[:2])

project = "FLUNED"
copyright = "2025, Marco De Pietri"
author = "Marco De Pietri"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output ------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
