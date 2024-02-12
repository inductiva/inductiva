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


def resolve_output_path(path: Optional[types.Path]) -> pathlib.Path:
    """Resolve a path relative to the output_dir

    Args:
        path: Path to a file or directory.
    """
    resolved_path = pathlib.Path.cwd()

    if inductiva.get_output_dir():
        resolved_path = pathlib.Path(inductiva.get_output_dir())

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


def _unzip(zip_path: pathlib.Path):
    """Unzip a zip archive and remove the .zip."""

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall()

    zip_path.unlink()


def download_from_url(url: str, unzip: bool = False) -> str:
    """Download a file from an URL.

    This function downloads from a set of files from an url.
    If the file is zipped, then the argument unzip allows to
    unzip everything exactly as it was zipped before.

    Args:
        url: The URL to download the file from.
        unzip: Whether to unzip the file after downloading.
    Returns:
        The path to the downloaded file.
    """
    # Get archive name from url passed
    local_path = url.split("/")[-1]
    local_path = pathlib.Path(local_path)

    try:
        logging.info("Downloading from URL to the local path: %s", local_path)
        downloaded_to, _ = urllib.request.urlretrieve(url, filename=local_path)
    except urllib.error.URLError as url_error:
        logging.error("Could not download file from %s", url)
        raise url_error

    # Unzip all files as they were zipped.
    if unzip and zipfile.is_zipfile(downloaded_to):
        local_path = local_path.with_suffix("")
        logging.info("Uncompressing the downloaded file to: %s", local_path)
        _unzip(downloaded_to)

    return str(local_path.absolute())
