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

def test_resolve_path():
    inductiva.working_dir = None
    resolved_path = files.resolve_path("protein.pdb")
    assert resolved_path == pathlib.Path.cwd().joinpath("protein.pdb")

    inductiva.working_dir = None
    resolved_path = files.resolve_path(None)
    assert resolved_path == pathlib.Path.cwd()

    inductiva.working_dir = "/tmp"
    resolved_path = files.resolve_path("protein.pdb")
    assert resolved_path == pathlib.Path("/tmp/protein.pdb")

    inductiva.working_dir = "/tmp"
    resolved_path = files.resolve_path("/protein.pdb")
    assert  resolved_path == pathlib.Path("/protein.pdb")

    inductiva.working_dir = "tmp"
    resolved_path = files.resolve_path("protein.pdb")
    assert resolved_path == pathlib.Path("tmp/protein.pdb")




