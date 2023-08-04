"""OpenFOAM module of the API for fluid dynamics."""
from typing import Optional, List
from uuid import UUID

from inductiva import types, tasks
from inductiva.simulation import Simulator


class OpenFOAM(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self, api_method: str = "fvm"):
        super().__init__()
        self.api_method = api_method + ".openfoam.run_simulation"

    @property
    def api_method_name(self) -> str:
        return self.api_method

    def run(
        self,
        input_dir: types.Path,
        commands: List[dict],
        resource_pool_id: Optional[UUID] = None,
        n_cores: int = 1,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            commands: List of commands to run using the OpenFOAM simulator.
            n_cores: Number of MPI cores to use for the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           resource_pool_id=resource_pool_id,
                           n_cores=n_cores,
                           commands=commands,
                           run_async=run_async)
