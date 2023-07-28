"""OpenFOAM module of the API for fluid dynamics."""
import pathlib
from typing import Optional, List, Union
from uuid import UUID

from inductiva import types
from inductiva.simulation import Simulator
from inductiva.tasks import Task


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
        run_async: bool = False,
    ) -> Union[pathlib.Path, Task]:
        """Run the simulation.

        Args:
            commands: List of commands to run using the OpenFOAM simulator.
            n_cores: Number of MPI cores to use for the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           output_dir=output_dir,
                           resource_pool_id=resource_pool_id,
                           n_cores=n_cores,
                           commands=commands,
                           run_async=run_async)
