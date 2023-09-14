"""Functions to download and load examples from URLs.

It contains several samples of molecules and a test_vehicle."""
import os
from typing import Optional
from urllib import request, error

from absl import logging

import inductiva


def load_insulin(save_dir: Optional[str] = None) -> str:
    """Load the insulin molecule.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1ZNI", save_dir)


def load_hemoglobin(save_dir: Optional[str] = None) -> str:
    """Load the hemoglobin molecule.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1A3N", save_dir)


def load_lysozyme(save_dir: Optional[str] = None) -> str:
    """Load lysozyme protein.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1AKI", save_dir)


def load_test_vehicle(save_dir: Optional[str] = None) -> str:
    """Load a test vehicle from Inductiva Github.

    Returns:
        The path to the downloaded OBJ file."""

    return download_inductiva_resources("test_vehicle.obj", save_dir)


def download_from_rcsb(pdb_id: str, save_dir: Optional[str] = None) -> str:
    """Download a PDB file from the RCSB database according 
    to its pdb ID.
    Args:
        pdb_id: The PDB identifier of the molecule to download.
    Returns:
        The path to the downloaded PDB file.
    """
    file_path = f"{pdb_id}.pdb"
    api_url = f"https://files.rcsb.org/download/{file_path}"

    file_path = download_from_external_url(api_url, file_path, save_dir)

    return file_path


def download_inductiva_resources(file_name: str, save_dir: Optional[str] = None) -> str:
    """Download a file from the resources folder of Inductiva Github.

    Args:
        file_name: The name of the file to download.
            Currently available: ["test_vehicle.obj", "alanine.pdb"]
    Returns:
        The path to the downloaded file.
    """

    api_url = "https://github.com/inductiva/inductiva/tree/" \
              f"main/resources/examples/{file_name}"
    file_path = download_from_external_url(api_url, file_name, save_dir)

    return file_path


def download_from_external_url(api_url: str, file_path: str, save_dir: Optional[str] = None) -> str:
    """Download a file from an external URL.

    Args:
        api_url: The URL to download the file from.
        file_path: The name of the file to download.
        save_dir: The directory to save the file to.
    Returns:
        The path to the downloaded file.
    """

    if save_dir is not None:
        save_dir = inductiva.utils.files.resolve_path(save_dir)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        file_path = os.path.join(save_dir, file_path)
    try:
        request.urlretrieve(api_url, file_path)
        logging.info("File downloaded to %s", file_path)
    except error.URLError as url_error:
        logging.error("Could not download file from %s", api_url)
        logging.error("%s", url_error)

    return  file_path
