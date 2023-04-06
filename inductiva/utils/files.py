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

    Args:
        path: Path to a file or directory. If it is relative, it is considered
            as being relative to the `inductiva.working_dir` directory. Else,
            the absolute path is returned. If `inductiva.working_dir` is None,
            then the current working directory is used.
        is_output_path: If True, the path is relative and
            `inductiva.working_dir` is not set, the path is considered as
            being relative to `inductiva.output_dir`.
    """
    root = pathlib.Path.cwd()

    if inductiva.working_dir is not None:
        root = pathlib.Path(inductiva.working_dir)
    else:
        if is_output_path:
            root = pathlib.Path(inductiva.output_dir)

    return pathlib.Path(root, path)
