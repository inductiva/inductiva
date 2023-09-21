"""Utilities for working with the file system."""
import os
import time
import pathlib
from typing import Optional
import urllib.error
import urllib.request
import zipfile

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


def _unzip(zip_path: pathlib.Path, dest_path: pathlib.Path):
    """Unzip a zip archive.

    If the ZIP contains one file, dest_path will be the path of that single
    file. If the ZIP contains multiple files, dest_path will be the path
    of the directory containing the files.
    """
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        files = zip_ref.namelist()
        if len(files) == 1:
            # Unzip the file to the directory of dest_path
            unzipped_file_path = pathlib.Path(
                zip_ref.extract(files[0], dest_path.parent))
            pathlib.Path(unzipped_file_path).replace(dest_path)
        else:
            dest_path.mkdir(parents=True, exist_ok=True)
            zip_ref.extractall(dest_path.parent)

    # Rename the file to dest_path
    zip_path.unlink()


def download_from_url(url: str, save_dir: Optional[str] = None) -> str:
    """Download a file from an URL.

    If the file is a ZIP archive containing a single file, the method
    will extract the file and rename it to `local_file_path`. If the ZIP
    contains more than one file, a ValueError will be raised.

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
    local_path = url.split("/")[-1]
    local_path = pathlib.Path(local_path)

    if save_dir is not None:
        local_path = pathlib.Path(save_dir, local_path)

    local_path = resolve_path(local_path)

    if not local_path.parent.exists():
        local_path.parent.mkdir(parents=True)

    try:
        downloaded_to, headers = urllib.request.urlretrieve(url)
    except urllib.error.URLError as url_error:
        logging.error("Could not download file from %s", url)
        raise url_error

    # File was downloaded to a temporary file
    downloaded_to = pathlib.Path(downloaded_to)
    is_zip = headers.get_content_type() == "application/zip"

    if is_zip:
        # Unzip the ZIP archive containing a single zip file to the
        # correct path

        if local_path.suffix == ".zip":  # Remove the .zip extension
            local_path = local_path.with_suffix("")

        _unzip(downloaded_to, local_path)
    else:
        # Rename the file to the correct path;
        # using replace instead of rename to have consistent behaviour
        # on windows and linux if the file already exists.
        # If the file exists, it will be overwritten.
        downloaded_to.replace(local_path)

    return str(local_path.absolute())
