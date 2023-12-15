import os

import pytest

from inductiva import mixins


@pytest.mark.parametrize("root_dir, expected_root_dir",
                         [("new_root_dir", "new_root_dir"),
                          ("new_root_dir", "new_root_dir")])
def test_set_root_dir(root_dir, expected_root_dir):
    file_manager = mixins.FileManager()
    file_manager.set_root_dir(root_dir)

    root_dir_path = file_manager.get_root_dir()
    name = root_dir_path.name.split("-")[0]
    assert name == str(expected_root_dir)


@pytest.mark.parametrize("source_file, render_args, target_file", [
    (os.path.join(os.path.dirname(__file__), "test_file.txt"), dict(), None),
    (os.path.join(os.path.dirname(__file__),
                  "test_file.txt"), dict(), "test2_file.txt"),
    (os.path.join(os.path.dirname(__file__),
                  "template_file.txt.jinja"), dict(text="yes"), None),
    (os.path.join(os.path.dirname(__file__), "template_file.txt.jinja"),
     dict(text="yes"), "test3_file.txt"),
])
def test_add_files(source_file, render_args, target_file):

    file_manager = mixins.FileManager()
    file_manager.set_root_dir(os.path.join(os.path.dirname(__file__), "assets"))
    target_file = file_manager.add_file(source_file, target_file, **render_args)

    assert os.path.isfile(os.path.join(file_manager.get_root_dir(),
                                       target_file))


@pytest.mark.parametrize("source_dir, render_args, target_dir", [
    (os.path.join(os.path.dirname(__file__), "test_dir"), dict(), None),
    (os.path.join(os.path.dirname(__file__), "test_dir"), dict(), "test2_dir"),
    (os.path.join(os.path.dirname(__file__),
                  "test_dir"), dict(text="yes"), None),
    (os.path.join(os.path.dirname(__file__),
                  "test_dir"), dict(text="yes"), "test3_dir"),
])
def test_add_dir(source_dir, render_args, target_dir):

    file_manager = mixins.FileManager()
    file_manager.set_root_dir(os.path.join(os.path.dirname(__file__), "assets"))

    target_dir = file_manager.add_dir(source_dir, target_dir, **render_args)
    assert os.path.isdir(os.path.join(file_manager.get_root_dir(),
                                      target_dir))
