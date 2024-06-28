"""Test the TemplateManager class."""
import tempfile
import pathlib
import shutil
import os

from jinja2 import exceptions
from pytest import mark
import pytest

from inductiva import TemplateManager
from inductiva.templating import manager

ASSETS_DIR = pathlib.Path(__file__).parent / "assets"


def _get_file_contents(file):
    with open(file, "r", encoding="utf-8") as f:
        return f.read()


@pytest.fixture(name="tmp_target_dir")
def fixture_tmp_target_dir():
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield pathlib.Path(tmpdirname)


@pytest.mark.parametrize("text", ["world", "inductiva"])
def test_render_dir__copies_and_renders_files(tmp_target_dir, text):
    # Determine if the directory structure is correctly copied and rendered
    # to the local directory.
    TemplateManager.render_dir(ASSETS_DIR,
                               tmp_target_dir,
                               overwrite=True,
                               text=text)
    assert os.path.isfile(tmp_target_dir / "template.txt")
    assert os.path.isfile(tmp_target_dir / "non_template.txt")
    assert os.path.isfile(tmp_target_dir / "folder/nested_template.txt")
    assert os.path.isfile(tmp_target_dir / "folder/nested_non_template.txt")

    expected_contents = f"hello {text}\n"

    assert _get_file_contents(tmp_target_dir /
                              "template.txt") == expected_contents
    assert _get_file_contents(tmp_target_dir /
                              "folder/nested_template.txt") == expected_contents


def test_render_dir__target_exists_exists__raises_exception(tmp_target_dir):
    # Determine if an exception is raised when any of the
    # target files exists and overwrite is set to False.

    # Copy a file to the target directory, in order to cause a conflict.
    shutil.copyfile(ASSETS_DIR / "non_template.txt",
                    tmp_target_dir / "non_template.txt")

    assert os.path.isfile(tmp_target_dir / "non_template.txt")

    with pytest.raises(FileExistsError):
        TemplateManager.render_dir(ASSETS_DIR,
                                   tmp_target_dir,
                                   overwrite=False,
                                   text="world")


def test_render_dir__target_dir_exists_overwites__renders_correctly(
        tmp_target_dir):
    # Determine if the directory is rendered when the target
    # exists and the overwrite flag is set.

    # Copy a file to the target directory, in order to cause a conflict.
    shutil.copyfile(ASSETS_DIR / "non_template.txt",
                    tmp_target_dir / "non_template.txt")
    assert os.path.isfile(tmp_target_dir / "non_template.txt")
    # With overwrite=True, the directory should be rendered.
    TemplateManager.render_dir(ASSETS_DIR,
                               tmp_target_dir,
                               overwrite=True,
                               text="world")


@mark.parametrize(
    "template_name, expected",
    [("template.txt.jinja", "hello world\n"),
     ("non_template.txt", "this is a non-template file\n"),
     ("folder/nested_template.txt.jinja", "hello world\n"),
     ("folder/nested_non_template.txt", "this is a non-template file\n")])
def test_render_file__correctly_renders_files(tmp_target_dir, template_name,
                                              expected):
    # Determine if template and non-template files are correctly rendered.

    TemplateManager.render_dir(ASSETS_DIR, tmp_target_dir, text="world")

    if manager.is_template(template_name):
        template_name = manager.strip_extension(template_name)
    current = _get_file_contents(tmp_target_dir / template_name)
    assert current == expected


@mark.parametrize("filename, is_template",
                  [("template.txt", False), ("template.txt.jinja", True),
                   ("folder/nested_template.txt", False),
                   ("folder/nested_template.txt.jinja", True)])
def test_is_template__correctly_identifies_templates(filename, is_template):
    # Determine if template files are correctly identified.
    assert manager.is_template(filename) == is_template


@mark.parametrize("filename, expected",
                  [("template.txt", "template.txt"),
                   ("template.txt.jinja", "template.txt"),
                   ("folder/template.txt", "folder/template.txt"),
                   ("folder/template.txt.jinja", "folder/template.txt")])
def test_strip_extension__correctly_removes_extension(filename, expected):
    # Determine if the template extension is correctly stripped from a
    # template filename.
    assert manager.strip_extension(filename) == expected


def test_render_file__missing_parameter__raises_exception(tmp_target_dir):
    # Determine if an exception is raised when a template is rendered but not
    # all required parameters are given.
    with pytest.raises(exceptions.UndefinedError):
        TemplateManager.render_dir(ASSETS_DIR, tmp_target_dir)
