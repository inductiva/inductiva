"""Utilities for working with the file system."""
import time
import pathlib
import inductiva

from inductiva import types


def get_timestamped_path(path: types.Path, sep: str = "-") -> pathlib.Path:
    """Return a path that does not exist by appending a timestamp.

    Args:
        path: Path to a file or directory.

    Returns:
        A path that does not exist by appending the timestamp.
    """
    path = pathlib.Path(path)
    timestamp = time.strftime("%Y-%m-%dT%Hh%Mm%Ss")

    name = f"{path.stem}{sep}{timestamp}"

    return path.with_name(name + path.suffix)


def resolve_path(path: types.Path) -> pathlib.Path:
    """Resolve a path relative to the Inductiva package working directory.

    Args:
        path: Path to a file or directory.
    """
    root = pathlib.Path.cwd()

    if inductiva.working_dir:
        root = pathlib.Path(inductiva.working_dir)

    return pathlib.Path(root, path)
