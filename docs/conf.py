# pylint: skip-file
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'Inductiva API Python client'
copyright = '2025, Inductiva Research Labs'
author = 'Inductiva Research Labs'

# Mock imports for modules that may not be available or cause issues
autodoc_mock_imports = [
    "setup",  # Prevents issues with `setup.py` executing sys.exit()
    "conftest",  # Avoids errors if pytest isn't installed
    "pytest",  # Mock pytest to avoid import failures in test-related files
]

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",  # Auto-generates docs from docstrings
    "sphinx.ext.napoleon",  # Supports Google/NumPy-style docstrings
    #    "sphinx.ext.viewcode",      # Adds links to source code
    "sphinx.ext.autosummary",  # Auto-generates a summary for modules
    'sphinx.ext.mathjax',
    'myst_parser',
    'sphinxcontrib.mermaid',
    'sphinx_copybutton',
    'sphinx_tabs.tabs',
    'sphinx_togglebutton',
    'sphinxcontrib.googleanalytics',
    'sphinxext.opengraph',
    'sphinx_sitemap',
    'sphinx_reredirects'
]

# Enable automatic docstring discovery
autosummary_generate = True

# Treat warnings as errors
nitpicky = True

autodoc_default_options = {
    "members": True,  # Include all public functions/methods
    "undoc-members": False,  # Include methods even if they lack docstrings?
    "private-members": False,  # Exclude private methods (_method_name)
    "special-members":
        "__init__",  # Ensure constructors (__init__) are documented
    "show-inheritance": True,  # Show class hierarchy
    "inherited-members": True,  # Show inherited methods from base classes
    "module-first": True,  # Show modules before class names
}

# Hides the full module name in class/method documentation
add_module_names = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
sphinx_tabs_valid_builders = ['html']

# The suffix(es) of source filenames.
# Note: important to list ipynb before md here: we have both md and ipynb
# copies of each notebook, and myst will choose which to convert based on
# the order in the source_suffix list. Notebooks which are not executed have
# outputs stored in ipynb but not in md, so we must convert the ipynb.
source_suffix = ['.rst', '.ipynb', '.md']

# The main toctree document.
main_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    'README.md',
    'markdown_sample.md',
    'task_state_diagram.md',
]

myst_links_external_new_tab = True

# Auto generate header anchors
myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'

html_theme_options = {
    "logo": {
        "image_light": "_static/inductiva-logo-black.svg",
        "image_dark": "_static/inductiva-logo-white.svg"
    },
    "collapse_navigation": False,
    "show_nav_level": 1
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_favicon = 'favicon.ico'

html_title = 'Python API for Cloud-based Simulation'

# Google Analytics
googleanalytics_id = "G-NHJ03C6M91"
googleanalytics_enabled = True

# OpenGraph protocol
ogp_site_name = "Inductiva.ai Docs"
ogp_site_url = "https://docs.inductiva.ai"
ogp_image = "https://docs.inductiva.ai/_static/inductiva-social-banner.jpg"

# sitemap.xml
# See https://sphinx-sitemap.readthedocs.io/
language = 'en'
version = 'local'
html_baseurl = 'https://docs.inductiva.ai/'

#For redirects to work there needs to be an actual html file
redirects = {
    "en/latest/how_to/run-parallel_simulations.html":
        "https://tutorials.inductiva.ai/how_to/run-parallel_simulations.html",
    "en/latest/how_to/manage_computational_resources.html":
        "https://tutorials.inductiva.ai/how_to/manage_computational_resources.html",
    "en/latest/how_to/set-up-elastic-machine-group.html":
        "https://tutorials.inductiva.ai/how_to/set-up-elastic-machine-group.html",
    "en/latest/how_to/set-up-mpi-cluster.html":
        "https://tutorials.inductiva.ai/how_to/set-up-elastic-machine-group.html",
    "en/latest/how_to/manage-remote-storage.html":
        "https://tutorials.inductiva.ai/how_to/set-up-elastic-machine-group.html",
    "en/latest/how_to/manage_tasks.html":
        "https://tutorials.inductiva.ai/how_to/set-up-elastic-machine-group.html",
    "en/latest/how_to/manage_and_retrieve_results.html":
        "https://tutorials.inductiva.ai/how_to/set-up-elastic-machine-group.html",
}
