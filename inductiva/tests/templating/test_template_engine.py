"""Test the TemplateManager class."""
import os
import pathlib
import shutil
import tempfile

import pytest

from inductiva import TemplateManager

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

    expected_contents = "hello " + text

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
