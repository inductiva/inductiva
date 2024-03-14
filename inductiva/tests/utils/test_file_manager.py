"""Tests for the file manager mixin."""
import pathlib
import os
import tempfile

import pytest
from pytest import mark

from inductiva.utils import file_manager

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
    suffix = file_manager._gen_unique_suffix(name, listnames)  #pylint: disable=protected-access
    suffix_tag = suffix.split(file_manager.SUFFIX_SEPARATOR)[-1]

    assert suffix_tag == expected_output


def test_set_root_dir__valid_input__creates_folder():
    """Check that with a non-empty input name the root dir is created.
    
    Goal: Test root dir creation from a given name. The folder is
    created in the current working dir, so we want to verify that the
    root dir was actually created, that the name is correct and that it
    is instantiated in the current working_dir.
    """

    root_dir = "test_root_dir"
    manager = file_manager.FileManager()
    manager.set_root_dir(root_dir)

    created_root_dir = manager.get_root_dir()
    assert os.path.isdir(created_root_dir)
    assert root_dir in created_root_dir.name
    assert pathlib.Path.cwd() == created_root_dir.parent


def test_set_root_dir_null_input__raises_error():
    """Check that with a null input name an error is raised.
    
    Goal: Test the usability of the set_root_dir method when root_dir is None.
    We capture the error that we raise (ValueError) and match the sentence to
    validate it is the correct error.
    """
    manager = file_manager.FileManager()
    with pytest.raises(ValueError) as excinfo:
        manager.set_root_dir(None)

    assert str(excinfo.value) == "Given root directory cannot be None"


def test_set_root_dir__root_exists__disable_filemanager_suffix(monkeypatch):
    """Test the env variable to disable the suffix creation.
    
    Goal: Test that setting the env variable to True, disables the auto-suffix
    and raises an error informing the user that the directory already exists.
    """
    monkeypatch.setenv("INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX", "True")
    root_dir = ASSETS_DIR

    with pytest.raises(FileExistsError) as excinfo:
        manager = file_manager.FileManager()
        manager.set_root_dir(root_dir)

    assert str(excinfo.value) == f"Directory {root_dir} already exists."


def test_get_root_dir__root_is_none__create_default_dir():
    """Check that with a null root dir the default dir is created.
    
    Goal: Test the scenario when the root_dir is created by default, when
    no root_dir is set first. Hence, the workflow should be still valid to
    not raise errors during the workflow. 
    """
    manager = file_manager.FileManager()
    original_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as temp_dir:
        os.chdir(temp_dir)
        root_dir = manager.get_root_dir()
        assert os.path.isdir(root_dir)
        assert root_dir.name == manager.DEFAULT_ROOT_DIR
        os.chdir(original_dir)


@mark.parametrize("source_file, target_file",
                  [(path("test_file.txt"), None),
                   (path("test_file.txt"), "copyed_file.txt")])
def test_copy_file__valid_file__to_target(source_file, target_file):
    """Check that a non-empty input file is copied to the given target.
    
    Goal: Test one of the functionalities of the FileManager().copy_file method,
    which copies a source file into a target file, without rendering any args.
    In this case, the copy_file function just copies the file directly.
    Hence, we verify if the target_file exists after being copied.
    """

    manager = file_manager.FileManager()
    target_file = manager.copy_file(source_file, target_file)
    print(target_file)

    root_dir = manager.get_root_dir()
    assert os.path.isfile(os.path.join(root_dir, target_file))


@mark.parametrize("source_dir, target_dir", [
    (ASSETS_DIR, None),
    (ASSETS_DIR, "copied_dir"),
])
def test_copy_dir__valid_dir__no_args(source_dir, target_dir):
    """Check that an existing input dir is copied to the given target dir.
    
    Goal: Test one of the functionalities of FileManager().copy_dir which just
    copies a source dir into a new target dir. If no arguments are passed, then
    no rendering occurs, which means that even template files are copied as is.
    All files are copied with the same exact name.
    """

    manager = file_manager.FileManager()
    target_dir = manager.copy_dir(source_dir, target_dir)

    root_dir = manager.get_root_dir()
    target_path = os.path.join(root_dir, target_dir)
    assert os.path.isdir(target_path)
    assert os.listdir(ASSETS_DIR) == os.listdir(target_path)


def test_validate_destination__raise_error__on_file_existence():
    """Check that the validation method raises an error when a file exists.
    
    Goal: Test if when another file with the same name already exists, the
    method raises a FileExistsError. This tests the overwrite mechanism of the
    FileManager.
    """

    target_file = os.path.join(ASSETS_DIR, "test_file.txt")
    manager = file_manager.FileManager()

    manager.set_root_dir("root_dir")
    manager.copy_file(target_file)

    with pytest.raises(FileExistsError) as excinfo:
        manager._check_precopy_dir(ASSETS_DIR, "root_dir")  # pylint: disable=protected-access

    assert "already exists" in str(excinfo.value)


def test_validate_destination__raise_error__on_dir_existence():
    """Check that the validation method raises an error when a file exists.
    
    Goal: Test if when a file with the same name of any file in the directory
    being passed already exists, the method raises a FileExistsError. This tests
    the overwrite mechanism of the FileManager.
    """

    manager = file_manager.FileManager()

    manager.set_root_dir("root_dir")
    manager.copy_dir(ASSETS_DIR)

    with pytest.raises(FileExistsError) as excinfo:
        manager._check_precopy_dir(ASSETS_DIR, "root_dir")  # pylint: disable=protected-access

    assert "already exists" in str(excinfo.value)
