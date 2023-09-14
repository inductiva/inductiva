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
        n_cores: int = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation. If None,
              the maximum number of cores available in the machine group will be
              used.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           input_filename=sim_config_filename,
                           n_cores=n_cores,
                           run_async=run_async)
