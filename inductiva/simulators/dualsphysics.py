"""DualSPHysics simulator module of the API."""

from typing import Optional, List

from inductiva import types, tasks, resources
from inductiva.simulators import Simulator


class DualSPHysics(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sph.dualsphysics.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: List[dict],
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[types.Path] = "",
    ) -> tasks.Task:
        """Executes a DualSPHysics simulation.

        Args:
            input_dir: Directory with simulation input files.
            sim_config_filename: Simulation config file.
            machine_group: Machine group for simulation.
            storage_dir: Directory for storing results.

        Returns:
            tasks.Task: An object representing the simulation task.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           commands=commands,
                           storage_dir=storage_dir)
