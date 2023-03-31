"""OpenFOAM module of the API for fluid dynamics."""
import pathlib
from typing import Optional

from inductiva import types
from inductiva.simulation import Simulator


class OpenFOAM(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "fvm.opemfoam.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        output_dir: Optional[types.Path] = None,
        openfoam_solver: str = "isoFoam",
        n_cores=1,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
            openfoam_solver: specific solver to simulate with OpenFOAM.
                OpenFOAM contains lots of solvers inside of it, which are used
                to call the run simulation through terminal, e.g.,
                [isoFoam, sonicFoam, ...]. The default solver is isoFoam.
        """
        return super().run(
            input_dir,
            output_dir,
            openfoam_solver=openfoam_solver,
            n_cores=n_cores,
        )
