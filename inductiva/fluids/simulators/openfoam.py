"""OpenFOAM module of the API for fluid dynamics."""
import pathlib
from typing import Optional

from inductiva.types import Path
from inductiva.fluids.simulators._simulator import Simulator


class OpenFOAM(Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "fvm.opemfoam.run_simulation"

    def simulate(self,
                 openfoam_solver: str,
                 output_dir: Optional[Path] = None,
                 n_cores=1) -> pathlib.Path:
        """Run the simulation.
        Args:
            n_cores: Number of MPI cores to use for the simulation.
            openfoam_solver: specific solver to simulate with OpenFOAM. 
            OpenFOAM contains lots of solvers inside of it, which are used
            to call the run simulation through terminal, e.g.,
            [isoFoam, sonicFoam,...] 
        """
        return super().simulate(output_dir=output_dir,
                                openfoam_solver=openfoam_solver,
                                n_cores=n_cores)
