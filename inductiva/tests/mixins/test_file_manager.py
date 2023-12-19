"""Tests for the file manager mixin."""
import pathlib
import io
import os

import pytest
from pytest import mark

from inductiva import mixins

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "assets")


def path(filename):
    """Get the path to a file in the assets directory."""
    return os.path.join(ASSETS_DIR, filename)


@mark.parametrize("name, listnames, expected_output",
                  [("test", ["inductiva", "."], ""),
                   ("test", ["inductiva", "test"], "2"),
                   ("test", ["inductiva", "test__2"], "3"),
                   ("test", ["inductiva", "test-2"], ""),
                   ("inductiva", ["is", "awesome"], ""),
                   ("inductiva", ["inductiva", "test"], "2")])
def test_gen_suffix__input_name__name_with_suffix(name, listnames,
                                                  expected_output):
    """Verifies if an input name is returned with a suffix.
    
    Goal: Test creation of suffixes when a name already exists in a list
    of names, with the following structure {name}#{number}.
    This suffix should increase with each repetition and with the suffix
    already present in the list.
    """
    suffix = mixins.file_manager._gen_unique_suffix(name, listnames)  #pylint: disable=protected-access
    suffix_tag = suffix.split(mixins.file_manager.SUFFIX_SEPARATOR)[-1]

    assert suffix_tag == expected_output


@mark.parametrize("source_file, render_args, expected_output",
                  [(path("template_file.txt.jinja"), {
                      "text": "yes"
                  }, "This is a test file, which yes is a jinja file!"),
                   (path("template.yaml.jinja"), {
                       "simulator": "openfoam"
                   }, "simulator: openfoam")])
def test_render_file__valid_file__rendered_file(source_file, render_args,
                                                expected_output):
    """Verify that an input file is rendered correctly.
    
    Goal: Test the rendering of a file with given render args. Starting from a
    given template file we render it from a set of args. The target file is
    tested to verify if it was rendered correctly from the args.
    """

    with io.StringIO() as f:
        mixins.file_manager.render_file(source_file, f, **render_args)
        f.seek(0)
        file_content = f.read()
        assert expected_output in file_content


def test_set_root_dir__valid_input__creates_folder():
    """Check that with a non-empty input name the root dir is created.
    
    Goal: Test root dir creation from a given name. The folder is
    created in the current working dir, so we want to verify that the
    root dir was actually created, that the name is correct and that it
    is instantiated in the current working_dir.
    """

    root_dir = "test_root_dir"
    file_manager = mixins.FileManager()

    file_manager.set_root_dir(root_dir)
    created_root_dir = file_manager.get_root_dir()
    assert os.path.isdir(created_root_dir)
    assert root_dir in created_root_dir.name
    assert pathlib.Path(os.getcwd()) == created_root_dir.parent


def test_set_root_dir_null_input__raises_error():
    """Check that with a null input name an error is raised.
    
    Goal: Test the usability of the set_root_dir method when root_dir is None.
    We capture the error that we raise (ValueError) and match the sentence to
    validate it is the correct error.
    """
    file_manager = mixins.FileManager()
    with pytest.raises(ValueError) as excinfo:
        file_manager.set_root_dir(None)

    assert str(excinfo.value) == "Given root directory cannot be None"


def test_set_root_dir__root_exists__disable_filemanager_suffix(monkeypatch):
    """Test the env variable to disable the suffix creation.
    
    Goal: Test that setting the env variable to True, disables the auto-suffix
    and raises an error informing the user that the directory already exists.
    """
    monkeypatch.setenv("INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX", "True")
    root_dir = ASSETS_DIR

    with pytest.raises(FileExistsError) as excinfo:
        file_manager = mixins.FileManager()
        file_manager.set_root_dir(root_dir)

    assert str(excinfo.value) == f"Directory {root_dir} already exists."


def test_get_root_dir__root_is_none__create_default_dir():
    """Check that with a null root dir the default dir is created.
    
    Goal: Test the scenario when the root_dir is created by default, when
    no root_dir is set first. Hence, the workflow should be still valid to
    not raise errors during the workflow. 
    """
    file_manager = mixins.FileManager()
    root_dir = file_manager.get_root_dir()

    assert os.path.isdir(root_dir)
    assert root_dir.name == file_manager.DEFAULT_ROOT_DIR


@mark.parametrize("source_file, target_file",
                  [(path("test_file.txt"), None),
                   (path("test_file.txt"), "added_file.txt")])
def test_add_file__valid_file__to_target(source_file, target_file):
    """Check that a non-empty input file is added to the given target.
    
    Goal: Test on of the functionalities of the FileManager().add_file method,
    which adds a source file into a target file, without rendering any args.
    In this case, the add_file function just copies the file directly.
    Hence, we verify if the target_file exists after being copied.
    """

    file_manager = mixins.FileManager()
    target_file = file_manager.add_file(source_file, target_file)

    root_dir = file_manager.get_root_dir()
    assert os.path.isfile(os.path.join(root_dir, target_file))


@mark.parametrize("source_file, render_args, target_file",
                  [(path("template_file.txt.jinja"), {
                      "text": "yes"
                  }, None),
                   (path("template_file.txt.jinja"), {
                       "text": "yes"
                   }, "rendered_test_file.txt")])
def test_add_files__valid_file__render_totarget(source_file, render_args,
                                                target_file):
    """Check that a non-empty input file is render to the given target.
    
    Goal: Test another functionality of the FileManager().add_file method, which
    renders template files from render args into a new target file. The render
    process automatically copies and renders the file on the spot. In case,
    target_file is None, then the target_file name will be the source file
    without the template extension.
    """

    file_manager = mixins.FileManager()
    target_file = file_manager.add_file(source_file, target_file, **render_args)
    root_dir = file_manager.get_root_dir()
    assert os.path.isfile(os.path.join(root_dir, target_file))


@mark.parametrize("source_dir, target_dir", [
    (ASSETS_DIR, None),
    (ASSETS_DIR, "added_dir"),
])
def test_add_dir__valid_dir__no_args(source_dir, target_dir):
    """Check that an existing input dir is added to the given target dir.
    
    Goal: Test one of the functionalities of FileManager().add_dir which just
    copies a source dir into a new target dir. If no arguments are passed, then
    no rendering occurs, which means that even template files are copied as is.
    All files are copied with the same exact name.
    """

    file_manager = mixins.FileManager()
    target_dir = file_manager.add_dir(source_dir, target_dir)

    root_dir = file_manager.get_root_dir()
    target_path = os.path.join(root_dir, target_dir)
    assert os.path.isdir(target_path)
    assert os.listdir(ASSETS_DIR) == os.listdir(target_path)


@mark.parametrize("source_dir, render_args, target_dir", [(ASSETS_DIR, {
    "text": "yes",
    "simulator": "openfoam"
}, None), (ASSETS_DIR, {
    "text": "yes",
    "simulator": "openfoam"
}, "rendered_dir")])
def test_add_dir__valid_dir__render_args(source_dir, render_args, target_dir):
    """Check that an input dir with template files is render to a target dir.
    
    Goal: Test another of the functionalities of FileManager().add_dir which
    copies a directory into another one and renders the template files from the
    given render_args. All files are copied, even the ones that aren't template
    files. Then, the template files are rendered and keep the same name but
    without the template extension.
    Hence, the rendering occurs only over files with the template extension and
    if render_args are passed. Further, if render_args doesn't contain all
    arguments that are on the template files, then those extra arguments are 
    rendered with None. We leave this to the responsability of the user.
    """
    file_manager = mixins.FileManager()
    target_dir = file_manager.add_dir(source_dir, target_dir, **render_args)

    root_dir = file_manager.get_root_dir()
    target_path = os.path.join(root_dir, target_dir)

    expected_render_files = [
        file.split(mixins.file_manager.TEMPLATE_EXTENSION)[0]
        for file in os.listdir(ASSETS_DIR)
    ]
    actual_render_files = os.listdir(target_path)
    assert len(actual_render_files) > 0
    assert sorted(expected_render_files) == sorted(actual_render_files)
