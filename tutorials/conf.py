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
copyright = '2024, Inductiva Research Labs'
author = 'Inductiva Research Labs'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.mathjax', 'sphinx.ext.napoleon',
    'myst_parser', 'sphinxcontrib.mermaid', 'sphinx_copybutton',
    'sphinx_tabs.tabs', 'sphinx_togglebutton', 'sphinxcontrib.googleanalytics',
    'sphinxext.opengraph', 'sphinx_sitemap', 'sphinx_reredirects'
]

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    # other MyST extensions
]

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
googleanalytics_id = "G-N6GVDZ3621"
googleanalytics_enabled = True

# OpenGraph protocol
ogp_site_name = "Inductiva.ai Tutorials"
ogp_site_url = "https://tutorials.inductiva.ai"
ogp_image = "https://tutorials.inductiva.ai/_static/inductiva-social-banner.jpg"

# sitemap.xml
# See https://sphinx-sitemap.readthedocs.io/
language = 'en'
version = 'local'
html_baseurl = 'https://tutorials.inductiva.ai/'

redirects = {
    "_sources/simulators/QuantumEspresso.md":
        "/simulators/QuantumEspresso.html",
    "_sources/simulators/Reef3D.md":
        "/simulators/Reef3D.html",
    "_sources/simulators/SCHISM.md":
        "/simulators/SCHISM.html",
    "_sources/simulators/XBeach.md":
        "/simulators/XBeach.html",
    "_sources/simulators/AmrWind.md":
        "/simulators/AmrWind.html",
    "_sources/simulators/OpenFOAM.md":
        "/simulators/OpenFOAM.html",
    "_sources/simulators/SPlisHSPlasH.md":
        "/simulators/SPlisHSPlasH.html",
    "_sources/simulators/FDS.md":
        "/simulators/FDS.html",
    "_sources/simulators/SWAN.md":
        "/simulators/SWAN.html",
    "_sources/simulators/CaNS.md":
        "/simulators/CaNS.html",
    "_sources/simulators/FVCOM.md":
        "/simulators/FVCOM.html",
    "_sources/simulators/SWASH.md":
        "/simulators/SWASH.html",
    "_sources/simulators/NWChem.md":
        "/simulators/NWChem.html",
    "_sources/simulators/GROMACS.md":
        "/simulators/GROMACS.html",
    "_sources/simulators/DualSPHysics.md":
        "/simulators/DualSPHysics.html",
    "_sources/simulators/overview.md":
        "/simulators/overview.html",
    "_sources/simulators/OpenFAST.md":
        "/simulators/OpenFAST.html",
    "_sources/intro_to_api/custom_docker_images.md":
        "/intro_to_api/custom_docker_images.html",
    "_sources/intro_to_api/computational-infrastructure.md":
        "/intro_to_api/computational-infrastructure.html",
    "_sources/intro_to_api/tasks.md":
        "/intro_to_api/tasks.html",
    "_sources/intro_to_api/data_flow.md":
        "/intro_to_api/data_flow.html",
    "_sources/intro_to_api/how_it_works.md":
        "/intro_to_api/how_it_works.html",
    "_sources/intro_to_api/configuring-simulators.md":
        "/intro_to_api/configuring-simulators.html",
    "_sources/pdes/heat-2-finite-differences.md":
        "/pdes/heat-2-finite-differences.html",
    "_sources/intro_to_api/projects.md":
        "/intro_to_api/projects.html",
    "_sources/pdes/heat-1-an-introduction.md":
        "/pdes/heat-1-an-introduction.html",
    "_sources/pdes/heat-3-PINN.md":
        "/pdes/heat-3-PINN.html",
    "_sources/pdes/heat-4-neurosolver.md":
        "/pdes/heat-4-neurosolver.html",
    "_sources/intro_to_api/shared_dedicated_resources.md":
        "/intro_to_api/shared_dedicated_resources.html",
    "_sources/intro_to_api/templating.md":
        "/intro_to_api/templating.html",
    "_sources/generating-synthetic-data/synthetic-data-generation-3.md":
        "/generating-synthetic-data/synthetic-data-generation-3.html",
    "_sources/generating-synthetic-data/synthetic-data-generation-1.md":
        "/generating-synthetic-data/synthetic-data-generation-1.html",
    "_sources/generating-synthetic-data/synthetic-data-generation-5.md":
        "/generating-synthetic-data/synthetic-data-generation-5.html",
    "_sources/generating-synthetic-data/synthetic-data-generation-4.md":
        "/generating-synthetic-data/synthetic-data-generation-4.html",
    "_sources/generating-synthetic-data/synthetic-data-generation-2.md":
        "/generating-synthetic-data/synthetic-data-generation-2.html",
    "_sources/index.md":
        "/index.html",
    "intro/_to/_api/templating.html":
        "/intro_to_api/templating.html",
    "en/latest/generating-synthetic-data/synthetic-data-generation-3.html":
        "/generating-synthetic-data/synthetic-data-generation-3.html",
    "en/latest/generating-synthetic-data/synthetic-data-generation-6.html":
        "/generating-synthetic-data/synthetic-data-generation-6.html",
    "en/latest/generating-synthetic-data/synthetic-data-generation-5.html":
        "/generating-synthetic-data/synthetic-data-generation-5.html",
    "en/latest/index.html":
        "/index.html",
    "en/latest/generating-synthetic-data/synthetic-data-generation-4.html":
        "/generating-synthetic-data/synthetic-data-generation-4.html",
    "en/latest/search.html":
        "/search.html",
    "en/latest/generating-synthetic-data/synthetic-data-generation-1.html":
        "/generating-synthetic-data/synthetic-data-generation-1.html",
    "en/latest/genindex.html":
        "/genindex.html",
    "en/latest/generating-synthetic-data/synthetic-data-generation-2.html":
        "/generating-synthetic-data/synthetic-data-generation-2.html",
}
