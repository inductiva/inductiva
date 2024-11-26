"""Utilities for working with the file system."""
import os
import time
import pathlib
from typing import Optional
import urllib.error
import urllib.request
import zipfile

import logging

import requests
from tqdm import tqdm

import inductiva


def get_timestamped_path(path_str: str, sep: str = "-") -> pathlib.Path:
    """Return a path that does not exist by appending a timestamp.

    Args:
        path: Path to a file or directory.

    Returns:
        A path that does not exist by appending the timestamp.
    """
    path = pathlib.Path(path_str)
    timestamp = time.strftime("%Y-%m-%dT%Hh%Mm%Ss")

    name = f"{path.stem}{sep}{timestamp}"

    return path.with_name(name + path.suffix)


def get_path_size(path_str: str) -> float:
    """Return the size of a path in bytes.
    
    Args:
        path: Path to a file or directory.
    
    Returns:
        The size of the path in bytes.
    """

    path = pathlib.Path(path_str)

    if not path.exists():
        raise FileNotFoundError(f"Path '{path}' does not exist.")

    if path.is_file():
        return path.stat().st_size

    size = 0

    for root, _, files in os.walk(path):
        root = pathlib.Path(root)
        size += root.stat().st_size
        for file in files:
            fp = root / file
            size += fp.stat().st_size
    return size


def resolve_output_path(path: Optional[str]) -> pathlib.Path:
    """Resolve a path relative to I/O directory

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


def _unzip(zip_path: pathlib.Path, unzip_path: Optional[pathlib.Path] = None):
    """Unzip a zip archive and remove the .zip."""

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        for member in tqdm(zip_ref.infolist(), desc="Extracting "):
            try:
                zip_ref.extract(member, path=unzip_path)
            except zipfile.error as e:
                print(f"Error extracting {member.filename}: {e}")
                pass

    zip_path.unlink()


def tqdm_refreshhook(tqdm_bar):

    def update_to(blocks, block_size, total_size):
        if total_size is not None:
            tqdm_bar.total = total_size
        tqdm_bar.n = blocks * block_size
        tqdm_bar.refresh()

    return update_to


def download_from_url(
    url: str,
    unzip: bool = False,
    path: Optional[str] = None,
) -> str:
    """Download a file from an URL.

    This function downloads from a set of files from an url.
    If the file is zipped, then the argument unzip allows to
    unzip everything exactly as it was zipped before.

    Args:
        url: The URL to download the file from.
        unzip: Whether to unzip the file after downloading.
        path: The directory path where the file should be saved. If None, 
            the file will be saved in the current working directory.

    Returns:
        The path to the downloaded file.
    """
    try:
        #Get head to check the filename
        filename = url.split("/")[-1]
        response = requests.head(url, allow_redirects=True, timeout=3)

        # Check if the 'Content-Disposition' header is present
        content_disposition = response.headers.get("Content-Disposition")
        if content_disposition and "filename=" in content_disposition:
            filename = content_disposition.split("filename=")[1].strip('"')

        logging.info("■ Downloading from URL to %s", filename)
        with tqdm(unit="B",
                  unit_scale=True,
                  unit_divisor=1024,
                  miniters=1,
                  desc=filename) as t:
            downloaded_to, _ = urllib.request.urlretrieve(
                url,
                filename=pathlib.Path(filename),
                reporthook=tqdm_refreshhook(t))
    except urllib.error.URLError as url_error:
        logging.error("Could not download file from %s", url)
        raise url_error

    resulting_path = downloaded_to

    # Unzip all files as they were zipped.
    if unzip and zipfile.is_zipfile(downloaded_to):
        resulting_path = (pathlib.Path(path)
                          if path else downloaded_to.with_suffix(""))
        logging.info("■ Uncompressing the downloaded archive to %s",
                     resulting_path)
        _unzip(downloaded_to, resulting_path if path else None)
        logging.info("")

    return str(resulting_path.absolute())
