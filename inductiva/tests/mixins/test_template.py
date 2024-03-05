import tempfile
import pathlib
import os

import pytest

from inductiva.mixins.template import Template


ASSETS_DIR = pathlib.Path(__file__).parent / "assets"


@pytest.fixture
def templatr():
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        yield Template(ASSETS_DIR)

def test_render_file__missing_parameters__raises_exception(templatr):
    with pytest.raises(Exception):
        templatr.render_file("template.txt.jinja")

def test_render_file__default_name__uses_template_name(templatr):
    # Determine if the template file is rendered to the local directory
    # with the same name as the template file, but with the template
    # extension stripped.
    templatr.render_file("template.txt.jinja", text="world")
    assert os.path.isfile("template.txt")

def test_render_file__dest_file_exists__raises_exception(templatr):
    # Determine if an exception is raised when the destination file exists.
    templatr.render_file("template.txt.jinja", text="world")
    with pytest.raises(FileExistsError):
        templatr.render_file("template.txt.jinja", text="world")

def test_render_file__dest_file_exists_overwrites__renders_correctly(templatr):
    # Determine if the destination file is overwritten when the
    # overwrite flag is set.
    templatr.render_file("template.txt.jinja", text="world")
    templatr.render_file("template.txt.jinja", text="world", overwrite=True)
    assert os.path.isfile("template.txt")


def test_render_file__nested_template__uses_template_name(templatr):
    # Determine if the nested template file is rendered to the local directory
    # with the same file name as the template file, but with the template
    # extension stripped.
    templatr.render_file("folder/nested_template.txt.jinja", text="world")
    assert os.path.isfile("nested_template.txt")

def test_render_file__nested_template__renders_to_name(templatr):
    # Determine if the nested template file is rendered to the default directory
    # with the given name.
    templatr.render_file("folder/nested_template.txt.jinja", "renamed.txt",
                         text="world")
    assert os.path.isfile("renamed.txt")

def test_render_file__nested_template__renders_to_name_in_dest_dir(templatr):
    # Determine if the nested template file is rendered to the given directory
    # with the given name.
    templatr.render_file("folder/nested_template.txt.jinja", "renamed.txt",
                            dest_dir="rendered",
                            text="world")
    assert os.path.isfile("rendered/renamed.txt")



def test_render_dir__default_dir__copies_and_renders_dir(templatr):
    # Determine if the directory structure is correctly copied and rendered
    # to the local directory.
    templatr.render_dir(text="world")
    assert os.path.isdir("empty")
    assert os.path.isfile("template.txt")
    assert os.path.isfile("non_template.txt")
    assert os.path.isfile("folder/nested_template.txt")
    assert os.path.isfile("folder/nested_non_template.txt")

def test_render_dir__dest_exists_exists__raises_exception(templatr):
    # Determine if an exception is raised when any of the
    # destination files exists.
    templatr.render_dir(text="world")
    with pytest.raises(FileExistsError):
        templatr.render_dir(text="world")

def test_render_dir__dest_dir_exists_overwites__renders_correctly(templatr):
    # Determine if the directory is rendered when the destination
    # exists and the overwrite flag is set.
    templatr.render_dir(text="world")
    templatr.render_dir(overwrite=True, text="world")

def test_render_dir__non_default_dest__renders_correctly(templatr):
    # Determine if the directory is rendered to a non-default destination.
    templatr.render_dir("rendered", text="world")
    assert os.path.isdir("rendered/empty")
    assert os.path.isfile("rendered/template.txt")
    assert os.path.isfile("rendered/non_template.txt")
    assert os.path.isfile("rendered/folder/nested_template.txt")
    assert os.path.isfile("rendered/folder/nested_non_template.txt")

def test_render_dir__dest_dir_exists__raises_exception(templatr):
    # Determine if an Exception is raised when the destination
    # exists and the overwrite flag is not set.
    templatr.render_dir("rendered", text="world")
    with pytest.raises(FileExistsError):
        templatr.render_dir("rendered", text="world")
