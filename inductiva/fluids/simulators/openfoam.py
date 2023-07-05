"""OpenFOAM module of the API for fluid dynamics."""
import pathlib
from typing import Optional, List
from uuid import UUID

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
        commands: List[dict],
        output_dir: Optional[types.Path] = None,
        resource_pool_id: Optional[UUID] = None,
        n_cores: int = 1,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            commands: List of commands to run using the OpenFOAM simulator.
            n_cores: Number of MPI cores to use for the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            output_dir=output_dir,
            resource_pool_id=resource_pool_id,
            n_cores=n_cores,
            commands=commands,
        )

    def run_async(
        self,
        input_dir: types.Path,
        method_name: str,
        resource_pool_id: Optional[UUID] = None,
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
                                 resource_pool_id=resource_pool_id,
                                 method_name=method_name,
                                 n_cores=n_cores,
                                 user_flags=openfoam_flags)
