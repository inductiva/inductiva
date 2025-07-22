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
from sphinxawesome_theme.postprocess import Icons

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'Inductiva.AI'
copyright = '2025, Inductiva.AI'
author = 'Inductiva.AI'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc', 'myst_parser', 'sphinx.ext.mathjax',
    'sphinx.ext.napoleon', 'sphinxcontrib.mermaid', 'sphinx_tabs.tabs',
    'sphinx_togglebutton', 'sphinxcontrib.googleanalytics',
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
html_permalinks_icon = Icons.permalinks_icon
html_theme = 'sphinxawesome_theme'

html_theme_options = {
    'show_prev_next': True,
    'show_scrolltop': True,
    'show_breadcrumbs': True,
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
shared_static_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "_shared_static"))
html_static_path = ['_static', shared_static_path]

html_css_files = ['css/custom.css', 'css/banner.css']
pygments_style = "monokai"

html_title = 'SCHISM'

# Google Analytics
googleanalytics_id = os.getenv("GTAG_WEBSITE", "GTM-K343XQD7")
googleanalytics_enabled = True

# OpenGraph protocol
ogp_site_name = "Inductiva.AI SCHISM"
ogp_site_url = "https://inductiva.ai/guides/schism"
ogp_image = "https://inductiva.ai/builds/schism/_static/inductiva-social-banner.jpg"

# sitemap.xml
# See https://sphinx-sitemap.readthedocs.io/
language = 'en'
version = 'local'
html_baseurl = 'https://inductiva.ai/guides/schism'

#save into static a js with the env var with the GTM code for the corrent env
#prod or dev
env_js_path = os.path.join(os.path.dirname(__file__), '_static', 'env.js')
os.makedirs(os.path.dirname(env_js_path), exist_ok=True)
with open(env_js_path, 'w') as f:
    f.write(f'window.env = {{ GTAG_WEBSITE: "{googleanalytics_id}" }};\n')
html_js_files = [
    'env.js',
    'discord.js',
    'gtm_func.js',
]

sys.path.insert(0, shared_static_path)


def setup(app):
    from banner_directive import BannerDirective
    app.add_directive("banner", BannerDirective)
    from banner_small_directive import BannerSmallDirective
    app.add_directive("banner_small", BannerSmallDirective)
