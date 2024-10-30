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

redirects = {
    "api_reference/computational_resources/compute_janitor.html":
        "/en/latest/api_reference/computational_resources/compute_janitor.html",
    "api_reference/computational_resources/elasticgroup_class.html":
        "/en/latest/api_reference/computational_resources/elasticgroup_class.html",
    "api_reference/computational_resources/index.html":
        "/en/latest/api_reference/computational_resources/index.html",
    "api_reference/computational_resources/machinegroup_class.html":
        "/en/latest/api_reference/computational_resources/machinegroup_class.html",
    "api_reference/computational_resources/mpicluster_class.html":
        "/en/latest/api_reference/computational_resources/mpicluster_class.html",
    "api_reference/faq.html":
        "/en/latest/api_reference/faq.html",
    "api_reference/glossary.html":
        "/en/latest/api_reference/glossary.html",
    "api_reference/troubleshooting.html":
        "/en/latest/api_reference/troubleshooting.html",
    "api_reference/uninstall_inductiva.html":
        "/en/latest/api_reference/uninstall_inductiva.html",
    "api_reference/user_quotas.html":
        "/en/latest/api_reference/user_quotas.html",
    "en/latest/api_reference/user_quotas.html":
        "/en/latest/api_reference/tiers_and_quotas.html",
    "cli/access-storage.html":
        "/en/latest/cli/access-storage.html",
    "cli/cli-overview.html":
        "/en/latest/cli/cli-overview.html",
    "cli/managing-resources.html":
        "/en/latest/cli/managing-resources.html",
    "cli/overview.html":
        "/en/latest/cli/overview.html",
    "cli/tracking-tasks.html":
        "/en/latest/cli/tracking-tasks.html",
    "computational_resources/index.html":
        "/en/latest/api_reference/computational_resources/index.html",
    "computational_resources/machinegroup_class.html":
        "/en/latest/api_reference/computational_resources/machinegroup_class.html",
    "en/latest/_sources/api_reference/user_quotas.md":
        "/en/latest/api_reference/tiers_and_quotas.html",
    #this here need to be updated to point to docs at a latter date
    "en/latest/_sources/explore_api/shared_dedicated_resources.md":
        "https://tutorials.inductiva.ai/intro_to_api/shared_dedicated_resources.html",
    "en/latest/explore_api/shared_dedicated_resources.html":
        "https://tutorials.inductiva.ai/intro_to_api/shared_dedicated_resources.html",
    "explore_api/shared_dedicated_resources.html":
        "https://tutorials.inductiva.ai/intro_to_api/shared_dedicated_resources.html",
    "explore_api/computational-infrastructure.html":
        "https://tutorials.inductiva.ai/intro_to_api/computational-infrastructure.html",
    "en/latest/explore_api/computational-infrastructure.html":
        "https://tutorials.inductiva.ai/intro_to_api/computational-infrastructure.html",
    "en/latest/explore_api/configuring-simulators.html":
        "https://tutorials.inductiva.ai/intro_to_api/configuring-simulators.html",
    "en/latest/explore_api/data_flow.html":
        "https://tutorials.inductiva.ai/intro_to_api/data_flow.html",
    "en/latest/explore_api/how_it_works.html":
        "https://tutorials.inductiva.ai/intro_to_api/how_it_works.html",
    "en/latest/explore_api/tasks.html":
        "https://tutorials.inductiva.ai/intro_to_api/tasks.html",
    "en/latest/explore_api/templating.html":
        "https://tutorials.inductiva.ai/intro_to_api/templating.html",
    "explore_api/configuring-simulators.html":
        "https://tutorials.inductiva.ai/intro_to_api/configuring-simulators.html",
    "explore_api/data_flow.html":
        "https://tutorials.inductiva.ai/intro_to_api/data_flow.html",
    "explore_api/how_it_works.html":
        "https://tutorials.inductiva.ai/intro_to_api/how_it_works.html",
    "explore_api/tasks.html":
        "https://tutorials.inductiva.ai/intro_to_api/tasks.html",
    "explore_api/templating.html":
        "thttps://tutorials.inductiva.ai/intro_to_api/templating.html",
    "introduction/templating.html":
        "https://tutorials.inductiva.ai/intro_to_api/templating.html",
    #This here needs to be updated to tutorials at a latter date
    "how_to/index.html":
        "https://docs.inductiva.ai/en/latest/how_to/index.html",
    "how_to/manage-remote-storage.html":
        "https://docs.inductiva.ai/en/latest/how_to/manage-remote-storage.html",
    "how_to/manage_and_retrieve_results.html":
        "https://docs.inductiva.ai/en/latest/how_to/manage_and_retrieve_results.html",
    "how_to/manage_computational_resources.html":
        "https://docs.inductiva.ai/en/latest/how_to/manage_computational_resources.html",
    "how_to/manage_tasks.html":
        "https://docs.inductiva.ai/en/latest/how_to/manage_tasks.html",
    "how_to/run-parallel_simulations.html":
        "https://docs.inductiva.ai/en/latest/how_to/run-parallel_simulations.html",
    "how_to/set-up-elastic-machine-group.html":
        "https://docs.inductiva.ai/en/latest/how_to/set-up-elastic-machine-group.html",
    "how_to/set-up-mpi-cluster.html":
        "https://docs.inductiva.ai/en/latest/how_to/set-up-mpi-cluster.html",
    #simulators
    "en/latest/simulators/SPlisHSPlasH.html":
        "https://tutorials.inductiva.ai/simulators/SPlisHSPlasH.html",
    "en/latest/simulators/SWAN.html":
        "https://tutorials.inductiva.ai/simulators/SWAN.html",
    "simulators/AmrWind.html":
        "https://tutorials.inductiva.ai/simulators/AmrWind.html",
    "simulators/CaNS.html":
        "https://tutorials.inductiva.ai/simulators/CaNS.html",
    "simulators/DualSPHysics.html":
        "https://tutorials.inductiva.ai/simulators/DualSPHysics.html",
    "simulators/FDS.html":
        "https://tutorials.inductiva.ai/simulators/FDS.html",
    "simulators/GROMACS.html":
        "https://tutorials.inductiva.ai/simulators/GROMACS.html",
    "simulators/NWChem.html":
        "https://tutorials.inductiva.ai/simulators/NWChem.html",
    "simulators/OpenFAST.html":
        "https://tutorials.inductiva.ai/simulators/OpenFAST.html",
    "simulators/OpenFOAM.html":
        "https://tutorials.inductiva.ai/simulators/OpenFOAM.html",
    "simulators/overview.html":
        "https://tutorials.inductiva.ai/simulators/overview.html",
    "simulators/Reef3D.html":
        "https://tutorials.inductiva.ai/simulators/Reef3D.html",
    "simulators/SCHISM.html":
        "https://tutorials.inductiva.ai/simulators/SCHISM.html",
    "simulators/SPlisHSPlasH.html":
        "https://tutorials.inductiva.ai/simulators/SPlisHSPlasH.html",
    "simulators/SWAN.html":
        "https://tutorials.inductiva.ai/simulators/SWAN.html",
    "simulators/SWASH.html":
        "https://tutorials.inductiva.ai/simulators/SWASH.html",
    "simulators/XBeach.html":
        "https://tutorials.inductiva.ai/simulators/XBeach.html",
}
