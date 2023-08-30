"""SWASH module of the API."""
from typing import Optional

from inductiva.simulation import Simulator
from inductiva import types, tasks, resources


class SWASH(Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sw.swash.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        machine_group: Optional[resources.MachineGroup] = None,
        n_cores: int = 1,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           input_filename=sim_config_filename,
                           n_cores=n_cores,
                           run_async=run_async)
