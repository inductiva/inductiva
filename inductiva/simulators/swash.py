"""SWASH module of the API."""
from typing import Optional

from inductiva import types, tasks, resources, simulators


class SWASH(simulators.Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sw.swash.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[types.Path] = "",
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            machine_group: Optional machine group to run the simulation on.
            storage_dir: Directory for storing simulation results.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           input_filename=sim_config_filename,
                           storage_dir=storage_dir)
