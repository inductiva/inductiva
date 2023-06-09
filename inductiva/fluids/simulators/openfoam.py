"""OpenFOAM module of the API for fluid dynamics."""
import pathlib
from typing import Optional

from inductiva import types
from inductiva.simulation import Simulator


class OpenFOAM(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self, api_method: str = "fvm.openfoam.run_simulation"):
        super().__init__()
        self.api_method = api_method

    @property
    def api_method_name(self) -> str:
        return self.api_method

    def run(
        self,
        input_dir: types.Path,
        method_name: str,
        output_dir: Optional[types.Path] = None,
        n_cores: int = 1,
        track_logs: bool = False,
        **openfoam_flags: Optional[str],
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
            method_name: OpenFOAM method to run. This involves commands
                from pre-processing, to solvers and post-processing.
            openfoam_flags: Flags to pass to the openfoam method_name.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           output_dir=output_dir,
                           track_logs=track_logs,
                           method_name=method_name,
                           n_cores=n_cores,
                           user_flags=openfoam_flags)

    def run_async(
        self,
        input_dir: types.Path,
        method_name: str,
        n_cores: int = 1,
        **openfoam_flags: Optional[str],
    ) -> str:
        """Run the simulation asynchronously.
        
        Args:
            n_cores: Number of MPI cores to use for the simulation.
            method_name: OpenFOAM method to run. This involves commands
                from pre-processing, to solvers and post-processing.
            openfoam_flags: Flags to pass to the openfoam method_name.
            other arguments: See the documentation of the base class.
            """

        return super().run_async(input_dir,
                                 method_name=method_name,
                                 n_cores=n_cores,
                                 user_flags=openfoam_flags)
