"""Functions to download and load example molecules for GROMACS simulations 
from the RCSB PDB database."""
import urllib
from absl import logging


def load_insulin() -> str:
    """Load the insulin molecule.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1ZNI")


def load_hemoglobin() -> str:
    """Load the hemoglobin molecule.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1A3N")


def load_lysozyme() -> str:
    """Load lysozyme protein.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1AKI")


def download_from_rcsb(pdb_id: str) -> str:
    """Download a PDB file from the RCSB database according 
    to its pdb ID.
    Args:
        pdb_id: The PDB identifier of the molecule to download.
    Returns:
        The path to the downloaded PDB file.
    """
    api_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    local_file_path = f"{pdb_id}.pdb"

    try:
        urllib.request.urlretrieve(api_url, local_file_path)
        logging.info("File downloaded to %s", local_file_path)
    except urllib.error.URLError as e:
        logging.error("Could not download file from %s", api_url)
        logging.error("%s", e)
    return local_file_path
