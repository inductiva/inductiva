"""Utilities for working with the file system."""
import os
import time
import pathlib
from typing import Optional
import urllib.error
import urllib.request

from absl import logging

import inductiva
from inductiva import types


def find_path_to_package(package_dir: str):
    """Find path to package directory in inductiva library in local computer.

    Args:
        package_dir: Name of package directory.
    """
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", package_dir))


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


def resolve_path(path: Optional[types.Path]) -> pathlib.Path:
    """Resolve a path relative to the Inductiva package working directory.

    This takes in a path and returns it resolved relative to the
    "inductiva.working_dir" setting.
    Note:
    If the working directory is not set, the current working directory of
    the script is used.
    If the path is absolute, it will override the working directory and be
    returned as is.
    If the path is None, the working directory will be returned.
    Check the tests in "test_files_utils.py" for examples.


    Args:
        path: Path to a file or directory.
    """
    resolved_path = pathlib.Path.cwd()

    if inductiva.working_dir:
        resolved_path = pathlib.Path(inductiva.working_dir)

    if path is not None:
        resolved_path = pathlib.Path(resolved_path, path)

    return resolved_path


def get_sorted_files(data_dir: str,
                     file_format: str = "name",
                     split_token: str = "_"):
    """Returns list of files sorted according to [file_key].

    Order a set of .format files of the form
    ['name_1.format', 'name_2.format',...,'name_10.format',
    ...,'name_n.format'].

    The default sorting methods for list, list.sort()
    or sorted(list), order 'name_10.format' before 'name_2.format',
    which is not representative of the time series.

    In this function we sort according to the number..
    """

    if not os.path.exists(data_dir):
        raise IOError(f"Directory '{data_dir}' does not exist.")

    # Get a list of the files in the data directory.
    files = os.scandir(data_dir)

    # The files have file_format extension.
    files = [
        file for file in files if pathlib.Path(file.path).suffix == file_format
    ]

    # Sort the files to be read according to [file_key].
    def get_alphanum_key(file):
        file_name = pathlib.Path(file.path).stem
        file_name_splits = file_name.split(split_token)
        file_key = file_name_splits[-1]
        return int(file_key)

    files = sorted(files, key=get_alphanum_key)

    return files


def download_from_url(url: str,
                      local_file_path: Optional[str] = None,
                      save_dir: Optional[str] = None) -> str:
    """Download a file from an URL.

    Args:
        url: The URL to download the file from.
        local_file_path: Name attributed to the file to download.
            If None is passed, it attributes the name given to
            the file on the url.
        save_dir: The directory to save the file to. If None
            is passed, this will download to the current working
            directory. If the save_dir passed does not exist, this
            method will try to create it. 
    Returns:
        The path to the downloaded file.
    """
    if local_file_path is None:
        local_file_path = url.split("/")[-1]

    if save_dir is not None:
        local_file_path = pathlib.Path(save_dir, local_file_path)

    local_file_path = resolve_path(local_file_path)

    if not local_file_path.parent.exists():
        local_file_path.parent.mkdir(parents=True)

    try:
        urllib.request.urlretrieve(url, local_file_path)
        logging.info("File downloaded to %s", local_file_path)
        return str(local_file_path.absolute())
    except urllib.error.URLError as url_error:
        logging.error("Could not download file from %s", url)
        raise url_error
