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


def resolve_path(path: types.Path, is_output_path=False) -> pathlib.Path:
    """Resolve a path relative to the Inductiva package working directory.

    The base for a provided relative path has the following precedence:
     1. `inductiva.working_dir`
     2. `inductiva.output_dir`
     3. The current working directory
    Where a smaller number indicates a higher precedence.

    Note that `inductiva.output_dir` is only used if `is_output_path` is True.

    Args:
        path: Path to a file or directory.
        is_output_path: If true, consider `inductiva.output_dir`.
    """
    root = pathlib.Path.cwd()

    if is_output_path and inductiva.output_dir:
        root = pathlib.Path(inductiva.output_dir)

    if inductiva.working_dir:
        root = pathlib.Path(inductiva.working_dir)

    return pathlib.Path(root, path)
