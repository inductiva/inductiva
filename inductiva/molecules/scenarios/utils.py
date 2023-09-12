"""Tools to analyze molecular dynamics simulations."""
try:
    import MDAnalysis as mda
    from MDAnalysis import transformations
    from MDAnalysis.analysis import align
except ImportError:
    mda = None
    transformations = None
    align = None

from inductiva.utils import optional_deps


@optional_deps.needs_molecules_extra_deps
def unwrap_trajectory(topology, trajectory):
    """Unwrap visualization of the trajectory to deal with
    Periodic Boundary Conditions.
    Args:
        topology: Path to the topology file.
        trajectory: Path to the trajectory file."""
    universe = mda.Universe(topology, trajectory, guess_bonds=True)
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
