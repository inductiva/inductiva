"""Functions to download and load example molecules for GROMACS simulations 
from the RCSB PDB database."""
import urllib


def load_hemoglobin():
    """Load the hemoglobin molecule.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1A3N")


def load_lysozyme():
    """Load lysozyme protein.
    Returns:
        The path to the downloaded PDB file."""
    return download_from_rcsb("1AKI")


def download_from_rcsb(pdb_id):
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
        print(f"File downloaded to {local_file_path}")
    except urllib.error.URLError as e:
        print(f"Failed to download the file: {e}")
    return local_file_path
