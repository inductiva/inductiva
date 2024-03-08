"""Test files module."""
import os
import glob
import pathlib
import inductiva
from inductiva.utils import files

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "assets")


def test_get_timestamped_path_with_ext(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    path = tmp_path / "file.txt"
    path.touch()

    timestamped_path = files.get_timestamped_path(path)

    assert not timestamped_path.exists()
    assert timestamped_path.name.startswith("file-")
    assert timestamped_path.name.endswith(".txt")


def test_map_dir_files():
    """Check if the files are correctly mapped to subdirs."""

    assets_map = {
        ".": {
            "files": [
                "test_file.txt",
                "template.yaml.jinja",
                "template_file.txt.jinja",
            ]
        },
        "test": {
            "files": ["file"]
        }
    }
    files_map = files.map_dir_files(ASSETS_DIR)

    for key, contents in files_map.items():
        assert assets_map[key]["files"].sort() == contents["files"].sort()


def test_map_dir_files__with_jinja_ext():
    """Check if the files are correctly mapped to subdirs with extension tag."""

    assets_map = {
        ".": {
            "files": ["test_file.txt",],
            ".jinja": [
                "template.yaml.jinja",
                "template_file.txt.jinja",
            ]
        },
        "test": {
            "files": ["file"],
            ".jinja": []
        }
    }

    files_map = files.map_dir_files(ASSETS_DIR, ".jinja")

    assert assets_map == files_map


def test_get_timestamp_path(tmp_path: pathlib.Path):
    """Check if multiple dirs are created with different names."""
    n = 10

    path = tmp_path / "output"

    for _ in range(n):
        path = files.get_timestamped_path(path)
        path.mkdir()

    assert len(glob.glob(str(path.parent / "output-*"))) == n


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
