import glob
import pathlib
from .files import get_timestamped_path


def test_get_timestamped_path_with_ext(tmp_path: pathlib.Path):
    """Test get_timestamped_path function."""
    path = tmp_path / "file.txt"
    path.touch()

    timestamped_path = get_timestamped_path(path)

    assert not timestamped_path.exists()
    assert timestamped_path.name.startswith("file_")
    assert timestamped_path.name.endswith(".txt")


def test_get_timestamp_path(tmp_path: pathlib.Path):
    N = 10

    path = tmp_path / "output"

    for _ in range(N):
        path = get_timestamped_path(path)
        path.mkdir()

    assert len(glob.glob(str(path.parent / "output_*"))) == N
