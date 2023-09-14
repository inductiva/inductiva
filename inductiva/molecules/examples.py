"""Functions to download and load example molecules for GROMACS simulations 
from the RCSB PDB database."""
import requests


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
    response = requests.get(api_url, timeout=10)
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        local_file_path = f"{pdb_id}.pdb"
        # Save the response content (PDB file) to the local file
        with open(local_file_path, "wb") as pdb_file:
            pdb_file.write(response.content)

        print(f"File downloaded to {local_file_path}")
    else:
        print("Failed to download the file.")
    return local_file_path
