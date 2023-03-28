"""Utilities for working with the file system."""
import time
import pathlib


def get_timestamped_path(path):
    """Return a path that does not exist by appending a timestamp.

    Args:
        path: Path to a file or directory.

    Returns:
        A path that does not exist by appending the timestamp.
    """
    path = pathlib.Path(path)
    timestamp = time.strftime("%Y%m%d-%H%M%S")

    name = f"{path.stem}_{timestamp}"

    return path.with_stem(name)
