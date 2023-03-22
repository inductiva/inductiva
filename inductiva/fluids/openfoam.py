"""OpenFOAM module of the API."""
import os
import pathlib

import inductiva
from inductiva.types import Path


class OpenFOAM:
    """Class to invoke a generic OpenFoam simulation on the API.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        openfoam_solver: specific solver to simulate with OpenFOAM. 
            OpenFOAM contains lots of solvers inside of it, which are used
            to call the run simulation through terminal, e.g.,
            [isoFoam, sonicFoam,...] 
    """

    def __init__(
        self,
        sim_dir: Path,
        openfoam_solver: str,
    ):
        self.openfoam_solver= openfoam_solver
        self.sim_dir = pathlib.Path(sim_dir)

        if not os.path.isdir(sim_dir):
            raise ValueError("The provided path is not a directory.")

    def simulate(self, n_cores=1, output_dir=None) -> Path:
        """Run the simulation.

        Args:
            n_cores: Number of CPU cores.
            output_dir: Directory where the generated files will be stored.
        """
        return inductiva.openfoam.run_simulation(self.sim_dir,
                                                 self.openfoam_solver,
                                                 n_cores=n_cores,
                                                 output_dir=output_dir)
