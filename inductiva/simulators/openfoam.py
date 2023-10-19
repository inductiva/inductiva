"""OpenFOAM module of the API for fluid dynamics."""
from typing import Optional, List

from inductiva import types, tasks, resources, simulators


class OpenFOAM(simulators.Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API."""

    def __init__(self, version: str = "foundation"):
        super().__init__()
        self.api_method_name = f"fvm.openfoam_{version}.run_simulation"

    def run(self,
            input_dir: types.Path,
            commands: List[dict],
            machine_group: Optional[resources.MachineGroup] = None
           ) -> tasks.Task:
        """Run the simulation.

        Args:
            commands: List of commands to run using the OpenFOAM simulator.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           commands=commands)
