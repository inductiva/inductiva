from io import StringIO
import pathlib

from jinja2 import exceptions
from pytest import mark
import pytest

from inductiva.mixins.template import JinjaRendererAdapter


ASSETS_DIR = pathlib.Path(__file__).parent / "assets"


@mark.parametrize(
        "template_name, expected",
        [("template.txt.jinja", "hello world"),
         ("non_template.txt", "this is a non-template file"),
         ("folder/nested_template.txt.jinja", "hello world"),
         ("folder/nested_non_template.txt", "this is a non-template file")])
def test_jinja_adapter_render(template_name, expected):
    # Determine if template and non-template files are correctly rendered.
    renderer = JinjaRendererAdapter(ASSETS_DIR)
    fout = StringIO()
    renderer.render(template_name, fout, text="world")
    assert fout.getvalue() == expected

@mark.parametrize("filename, is_template",
                  [("template.txt", False),
                   ("template.txt.jinja", True),
                   ("folder/nested_template.txt", False),
                   ("folder/nested_template.txt.jinja", True)])
def test_jinja_adapter__is_adapter(filename, is_template):
    # Determine if the adapter correctly identifies template files.
    assert JinjaRendererAdapter.is_template(filename) == is_template

@mark.parametrize("filename, expected",
                  [("template.txt", "template.txt"),
                   ("template.txt.jinja", "template.txt"),
                   ("folder/template.txt", "folder/template.txt"),
                   ("folder/template.txt.jinja", "folder/template.txt")])
def test_jinja_adapter__strip_extension(filename, expected):
    # Determine if the adapter correctly strips the extension
    # from a template filename.
    assert JinjaRendererAdapter.strip_extension(filename) == expected

def test_jinja_adapter__missing_parameter__raises_exception():
    # Determine if an exception is raised when a template is rendered but not
    # all required parameters are given.
    renderer = JinjaRendererAdapter(ASSETS_DIR)
    with pytest.raises(exceptions.UndefinedError):
        renderer.render("template.txt.jinja", StringIO())

def test_jinja_adapter__invalid_template__raises_exception():
    # Determine if an exception is raised when a non-existing template is
    # rendered.
    renderer = JinjaRendererAdapter(ASSETS_DIR)
    with pytest.raises(exceptions.TemplateNotFound):
        renderer.render("non_existing_template", StringIO())