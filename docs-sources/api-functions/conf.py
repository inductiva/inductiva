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

# Customizations for sphinx-argparse-cli generated HTML file.
import re
import bs4
import pathlib

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'Inductiva.AI'
copyright = '2025, Inductiva.AI'
author = 'Inductiva.AI'

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
    'sphinx_reredirects',
    'sphinx_argparse_cli',
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

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    # other MyST extensions
]

# Hides the full module name in class/method documentation
add_module_names = False

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

html_css_files = ['css/custom.css', 'css/enable_sidebar_focus.css']
pygments_style = "monokai"

# SEO - Add any paths that contain templates here, relative to this directory.
templates_path = [
    os.path.relpath(
        os.path.join(os.path.dirname(__file__), "..",
                     "_shared_templates/_templates"))
]

html_title = "API Functions"
html_context = {
    "project_name":
        "Inductiva.AI",
    "project_description":
        "Inductiva.AI API Functions",
    "project_url":
        "https://inductiva.ai/guides/api-functions",
    "ogp_image":
        "https://inductiva.ai/builds/api-functions/_static/inductiva-social-banner.jpg",
    "keywords":
        "api-functions, Inductiva.AI"
}

# Google Analytics
googleanalytics_id = os.getenv("GTAG_WEBSITE")
googleanalytics_enabled = True

# sitemap.xml
# See https://sphinx-sitemap.readthedocs.io/
language = 'en'
version = 'local'
html_baseurl = 'https://inductiva.ai/guides/api-functions'

# Customizations for sphinx-argparse-cli generated HTML file.
# Acts as a plugin to extend the functionality of sphinx-argparse-cli.
def fix_newlines(text):
    return re.sub(r'(?<!\n)\n(?!\n)', ' ', text)

def fix_paragraphs(soup):
    for section in soup.select('section[id^="inductiva-"]'):
        for paragraph in section.find_all("p"):
            paragraph.string = fix_newlines(paragraph.get_text())

# Regex to match backtick-wrapped commands starting with `inductiva`
RE_BACKTICK_WRAPPED = re.compile(r"`(inductiva(?:\s+\w+)+)`")

def sub_backtick_wrapped(text):
    re_code = r"<code>\1</code>"
    return RE_BACKTICK_WRAPPED.sub(re_code, text)


def replace_backtick_wrapped(soup):
    for paragraph in soup.find_all("p"):
        if not paragraph.string:
            continue
        new_html = sub_backtick_wrapped(paragraph.string)
        paragraph.clear()
        paragraph.append(bs4.BeautifulSoup(new_html, "html.parser"))


def tidy_cli_html(file_path):
    text = file_path.read_text(encoding="utf-8")
    soup = bs4.BeautifulSoup(text, "html.parser")

    fix_paragraphs(soup)
    replace_backtick_wrapped(soup)

    file_path.write_text(str(soup), encoding="utf-8")

def on_build_finished(app, exception):
    output_dir = pathlib.Path(app.outdir)

    cli_html_dir = output_dir / "auto-cli"
    for html_file in cli_html_dir.glob("**/*.html"):
        tidy_cli_html(html_file)

def setup(app):
    app.connect("build-finished", on_build_finished)