"""Test files module."""
import json
import os
import glob
import pathlib
import zipfile

import inductiva
from inductiva.utils import files, data

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "assets")


def test_get_timestamped_path_with_ext(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    path = tmp_path / "file.txt"
    path.touch()

    timestamped_path = files.get_timestamped_path(str(path))

    assert not timestamped_path.exists()
    assert timestamped_path.name.startswith("file-")
    assert timestamped_path.name.endswith(".txt")


def test_get_timestamp_path(tmp_path: pathlib.Path):
    """Check if multiple dirs are created with different names."""
    # 10 was failing on Windows, due to path being too long
    n = 3

    path = tmp_path / "output"

    for _ in range(n):
        path = files.get_timestamped_path(str(path))
        path.mkdir()
    assert len(glob.glob(str(path.parent / "output-*"))) == n


def test_get_path_size__input_file(tmp_path: pathlib.Path):
    """Check if the size of a file is correctly calculated."""
    path = tmp_path / "file.txt"
    path.write_text("Hello, World!")

    size = files.get_path_size(str(path))

    assert pathlib.Path(path).stat().st_size == size


def test_get_path_size__input_dir(tmp_path: pathlib.Path):
    """Check if the size of a directory is correctly calculated."""
    expected_size = 0
    path_dir1 = tmp_path / "dir1"
    path_dir1.mkdir(parents=True, exist_ok=True)

    path_1 = path_dir1 / "file1.txt"
    (path_1).write_text("Hello, World!")

    path_2 = tmp_path / "file2.txt"
    (path_2).write_text("I'm not a very long file at all!")

    expected_size += tmp_path.stat().st_size
    expected_size += path_dir1.stat().st_size
    expected_size += path_1.stat().st_size
    expected_size += path_2.stat().st_size

    size = files.get_path_size(str(tmp_path))

    assert size == expected_size


def test_resolve_output_path():
    inductiva.set_output_dir(None)
    resolved_path = files.resolve_output_path("protein.pdb")
    assert resolved_path == pathlib.Path.cwd().joinpath("protein.pdb")

    resolved_path = files.resolve_output_path(None)
    assert resolved_path == pathlib.Path.cwd()

    inductiva.set_output_dir("/tmp")
    resolved_path = files.resolve_output_path("/protein.pdb")
    assert resolved_path == pathlib.Path("/protein.pdb")

    inductiva.set_output_dir("/tmp")
    resolved_path = files.resolve_output_path("protein.pdb")
    assert resolved_path == pathlib.Path("/tmp/protein.pdb")


def test_pack_input(tmp_path: pathlib.Path):
    """Test pack_input function."""
    dir_path = tmp_path / "dir"
    dir_path.mkdir()
    (dir_path / "file_in_dir.txt").write_text("Hello, Directory!")
    zip_name_dir = "test_zip_dir"

    zip_path_dir = data.pack_input(dir_path, zip_name_dir)
    assert zipfile.is_zipfile(zip_path_dir)

    with zipfile.ZipFile(zip_path_dir, "r") as zip_f:
        assert "sim_dir/file_in_dir.txt" in zip_f.namelist()


def test_extract_subdir_files(tmp_path: pathlib.Path):
    """Test extract_subdir_files function."""
    zip_path = tmp_path / "test.zip"
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    with zipfile.ZipFile(zip_path, "w") as zip_f:
        zip_f.writestr("artifacts/file.txt", "Hello, World!")

    with zipfile.ZipFile(zip_path, "r") as zip_f:
        data.extract_subdir_files(zip_f, "artifacts", output_dir)

    assert (output_dir / "file.txt").exists()


def test_zip_dir(tmp_path: pathlib.Path):
    """Test zip_dir function."""
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    (src_dir / "file.txt").write_text("Hello, World!")

    zip_name = "test_zip"
    zip_path = data.zip_dir(src_dir, zip_name)

    assert zipfile.is_zipfile(zip_path)


def test_uncompress_task_outputs(tmp_path: pathlib.Path):
    """Test uncompress_task_outputs function."""
    zip_path = tmp_path / "test.zip"
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    with zipfile.ZipFile(zip_path, "w") as zip_f:
        zip_f.writestr("artifacts/file.txt", "Hello, World!")

    data.decompress_zip(zip_path, output_dir)

    assert (output_dir / "artifacts/file.txt").exists()


def test_extract_zip_file_to_output_dir(tmp_path: pathlib.Path):
    """Test extract_zip_file_to_output_dir function."""
    zip_path = tmp_path / "test.zip"
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    files_list = ["file1.txt", "file2.txt", "file3.txt"]

    with zipfile.ZipFile(zip_path, "w") as zip_f:
        for file in files_list:
            zip_f.writestr(file, "Hello, World!")

    with zipfile.ZipFile(zip_path, "r") as zip_f:

        for file in files_list:
            #pylint: disable=protected-access
            data._extract_zip_file_to_dir(dest_dir=output_dir,
                                          remove_zip_file=zip_f,
                                          zip_path=file,
                                          filename=file)
            #pylint: enable=protected-access

    assert (output_dir / "file1.txt").exists()
    assert (output_dir / "file2.txt").exists()
    assert (output_dir / "file3.txt").exists()
