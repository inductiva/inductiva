"""OpenFOAM module of the API for fluid dynamics."""
import pathlib
from typing import Optional

from inductiva import types
from inductiva.simulation import Simulator


class OpenFOAM(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "fvm.openfoam.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        output_dir: Optional[types.Path] = None,
        openfoam_solver: str = "interFoam",
        n_cores: int = 1,
        track_logs: bool = False,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
            openfoam_solver: specific solver to simulate with OpenFOAM.
                OpenFOAM contains lots of solvers inside of it, which are used
                to call the run simulation through terminal, e.g.,
                [isoFoam, sonicFoam, ...]. The default solver is interFoam.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            output_dir=output_dir,
            track_logs=track_logs,
            openfoam_solver=openfoam_solver,
            n_cores=n_cores,
        )

    def run_async(
        self,
        input_dir: types.Path,
        openfoam_solver: str = "interFoam",
        n_cores: int = 1,
    ) -> str:
        """Run the simulation asynchronously.
        
        Args:
            n_cores: Number of MPI cores to use for the simulation.
            openfoam_solver: specific solver to simulate with OpenFOAM.
                OpenFOAM contains lots of solvers inside of it, which are used
                to call the run simulation through terminal, e.g.,
                [isoFoam, sonicFoam, ...]. The default solver is interFoam.
            other arguments: See the documentation of the base class.
            """

        return super().run_async(input_dir,
                                 openfoam_solver=openfoam_solver,
                                 n_cores=n_cores)
