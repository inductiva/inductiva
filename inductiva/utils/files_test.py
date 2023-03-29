"""Test files module."""
import glob
import pathlib
from .files import get_timestamped_path


def test_get_timestamped_path_with_ext(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    path = tmp_path / "file.txt"
    path.touch()

    timestamped_path = get_timestamped_path(path)

    assert not timestamped_path.exists()
    assert timestamped_path.name.startswith("file-")
    assert timestamped_path.name.endswith(".txt")


def test_get_timestamp_path(tmp_path: pathlib.Path):
    """Check if multiple dirs are created with different names."""
    n = 10

    path = tmp_path / "output"

    for _ in range(n):
        path = get_timestamped_path(path)
        path.mkdir()

    assert len(glob.glob(str(path.parent / "output-*"))) == n
