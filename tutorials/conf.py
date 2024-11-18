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
    "intro_to_api/how_it_works.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/how_it_works.html",
    "intro_to_api/tasks.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/tasks.html",
    "intro_to_api/shared_dedicated_resources.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/shared_dedicated_resources.html",
    "intro_to_api/data_flow.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/data_flow.html",
    "intro_to_api/computational-infrastructure.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/computational-infrastructure.html",
    "intro_to_api/computational-infrastructure.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/computational-infrastructure.html",
    "intro_to_api/templating.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/templating.html",
    "intro_to_api/configuring-simulators.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/configuring-simulators.html",
    "intro_to_api/projects.html":
        "https://docs.inductiva.ai/en/latest/intro_to_api/projects.html"
}
