"""OpenFOAM module of the API for fluid dynamics."""
from typing import Optional, List

from inductiva import types, tasks, resources, simulators


class OpenFOAM(simulators.Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "fvm.openfoam.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: List[dict],
        machine_group: Optional[resources.MachineGroup] = None,
        n_cores: int = 2,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            commands: List of commands to run using the OpenFOAM simulator.
            n_cores: Number of MPI cores to use for the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           n_cores=n_cores,
                           commands=commands,
                           run_async=run_async)
