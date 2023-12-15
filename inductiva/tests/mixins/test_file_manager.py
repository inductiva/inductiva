"""Tests for the file manager mixin."""
import os

import pytest

from inductiva import mixins


@pytest.mark.parametrize("root_dir", [("source_root_dir"), ("source_root_dir")])
def test_set_root_dir(root_dir):
    file_manager = mixins.FileManager()
    file_manager.set_root_dir(root_dir)

    root_dir_path = file_manager.get_root_dir()

    # Test creating multiple root directiories with the same name
    assert root_dir in root_dir_path.name


@pytest.mark.parametrize("source_file, render_args, target_file", [
    (os.path.join(os.path.dirname(__file__), "test_file.txt"), {}, None),
    (os.path.join(os.path.dirname(__file__),
                  "test_file.txt"), {}, "test2_file.txt"),
    (os.path.join(os.path.dirname(__file__), "template_file.txt.jinja"), {
        "text": "yes"
    }, None),
    (os.path.join(os.path.dirname(__file__), "template_file.txt.jinja"), {
        "text": "yes"
    }, "test3_file.txt"),
])
def test_add_files(source_file, render_args, target_file):

    file_manager = mixins.FileManager()
    file_manager.set_root_dir(os.path.join(os.path.dirname(__file__), "assets"))
    target_file = file_manager.add_file(source_file, target_file, **render_args)

    assert os.path.isfile(os.path.join(file_manager.get_root_dir(),
                                       target_file))


@pytest.mark.parametrize("source_dir, render_args, target_dir", [
    (os.path.join(os.path.dirname(__file__), "test_dir"), {}, None),
    (os.path.join(os.path.dirname(__file__), "test_dir"), {}, "test2_dir"),
    (os.path.join(os.path.dirname(__file__), "test_dir"), {
        "text": "yes"
    }, None),
    (os.path.join(os.path.dirname(__file__), "test_dir"), {
        "text": "yes"
    }, "test3_dir"),
])
def test_add_dir(source_dir, render_args, target_dir):

    file_manager = mixins.FileManager()
    file_manager.set_root_dir(os.path.join(os.path.dirname(__file__), "assets"))

    target_dir = file_manager.add_dir(source_dir, target_dir, **render_args)
    assert os.path.isdir(os.path.join(file_manager.get_root_dir(), target_dir))
