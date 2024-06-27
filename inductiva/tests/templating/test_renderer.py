"""Test for the renderers available."""
from io import StringIO
import pathlib

from jinja2 import exceptions
from pytest import mark
import pytest

from inductiva.templating.renderers import JinjaRenderer

ASSETS_DIR = pathlib.Path(__file__).parent / "assets"


@mark.parametrize(
    "template_name, expected",
    [("template.txt.jinja", "hello world"),
     ("non_template.txt", "this is a non-template file"),
     ("folder/nested_template.txt.jinja", "hello world"),
     ("folder/nested_non_template.txt", "this is a non-template file")])
def test_jinja_renderer_render(template_name, expected):
    # Determine if template and non-template files are correctly rendered.
    renderer = JinjaRenderer(ASSETS_DIR)
    fout = StringIO()
    renderer.render_file(template_name, fout, text="world")
    print(fout.getvalue())
    assert fout.getvalue() == expected


@mark.parametrize("filename, is_template",
                  [("template.txt", False), ("template.txt.jinja", True),
                   ("folder/nested_template.txt", False),
                   ("folder/nested_template.txt.jinja", True),
                   (pathlib.Path("folder/nested_template.txt.jinja"), True),
                   (pathlib.Path("folder/nested_template.txt"), False)])
def test_jinja_renderer__is_template(filename, is_template):
    # Determine if the renderer correctly identifies template files.
    renderer = JinjaRenderer(ASSETS_DIR)
    assert renderer.is_template(filename) == is_template


@mark.parametrize("filename, expected",
                  [("template.txt", "template.txt"),
                   ("template.txt.jinja", "template.txt"),
                   ("folder/template.txt", "folder/template.txt"),
                   ("folder/template.txt.jinja", "folder/template.txt")])
def test_jinja_renderer__strip_extension(filename, expected):
    # Determine if the renderer correctly strips the extension
    # from a template filename.
    file_dir = pathlib.Path(filename).parent
    renderer = JinjaRenderer(file_dir)
    assert renderer.strip_extension(filename) == expected


def test_jinja_adapter__missing_parameter__raises_exception():
    # Determine if an exception is raised when a template is rendered but not
    # all required parameters are given.
    file = "template.txt.jinja"
    renderer = JinjaRenderer(ASSETS_DIR)
    with pytest.raises(exceptions.UndefinedError):
        renderer.render_file(file, StringIO())


def test_jinja_renderer__invalid_template__raises_exception():
    # Determine if an exception is raised when a non-existing template is
    # rendered.
    renderer = JinjaRenderer("something")
    with pytest.raises(exceptions.TemplateNotFound):
        renderer.render_file("non_existing_template", StringIO())
