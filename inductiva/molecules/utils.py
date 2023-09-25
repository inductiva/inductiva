"""Tools to analyze molecular dynamics simulations."""
from typing import Optional

try:
    import MDAnalysis as mda
    from MDAnalysis import transformations
    from MDAnalysis.analysis import align
except ImportError:
    mda = None
    transformations = None
    align = None

from inductiva.utils import optional_deps, files


@optional_deps.needs_molecules_extra_deps
def unwrap_trajectory(topology, trajectory, guess_bonds=True):
    """Unwrap visualization of the trajectory to deal with
    Periodic Boundary Conditions.
    Args:
        topology: Path to the topology file.
        trajectory: Path to the trajectory file.
        guess_bonds: Whether to guess bonds between atoms."""
    universe = mda.Universe(topology, trajectory, guess_bonds=guess_bonds)
    atoms = universe.atoms
    transformation = transformations.unwrap(atoms)
    universe.trajectory.add_transformations(transformation)
    return universe


@optional_deps.needs_molecules_extra_deps
def align_trajectory_to_average(universe, trajectory_output_path):
    """Align the trajectory to the average structure.
    Args:
        universe: The universe MDAnalysis object.
        trajectory_output_path: Path to the aligned trajectory file."""
    average = align.AverageStructure(universe,
                                     universe,
                                     select="protein and name CA",
                                     ref_frame=0).run()
    average_trajectory = average.results.universe

    align.AlignTraj(universe,
                    average_trajectory,
                    select="protein and name CA",
                    filename=trajectory_output_path,
                    in_memory=False).run()


def download_pdb_from_rcsb(pdb_id: str, save_dir: Optional[str] = None) -> str:
    """Download a PDB file from the RCSB database according
    to its pdb ID.
    Args:
        pdb_id: The PDB identifier of the molecule to download.
        save_dir: The directory to save the file to. If None
            is passed, this will download to the current working
            directory. If the save_dir passed does not exist, this
            method will try to create it.
    Returns:
        The path to the downloaded PDB file.
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    file_path = files.download_from_url(pdb_url, save_dir)

    return file_path
