"""Test files module."""
import glob
import pathlib
import inductiva
from inductiva.utils import files


def test_get_timestamped_path_with_ext(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    path = tmp_path / "file.txt"
    path.touch()

    timestamped_path = files.get_timestamped_path(path)

    assert not timestamped_path.exists()
    assert timestamped_path.name.startswith("file-")
    assert timestamped_path.name.endswith(".txt")


def test_get_timestamp_path(tmp_path: pathlib.Path):
    """Check if multiple dirs are created with different names."""
    n = 10

    path = tmp_path / "output"

    for _ in range(n):
        path = files.get_timestamped_path(path)
        path.mkdir()

    assert len(glob.glob(str(path.parent / "output-*"))) == n


def test_get_path_size__input_file(tmp_path: pathlib.Path):
    """Check if the size of a file is correctly calculated."""
    path = tmp_path / "file.txt"
    path.write_text("Hello, World!")

    size = files.get_path_size(path)

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

    size = files.get_path_size(tmp_path)

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
