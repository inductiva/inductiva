"""FDS simulator module of the API."""

from typing import Optional

from inductiva import types, tasks, resources
from inductiva.simulators import Simulator


class FDS(Simulator):
    """Class to invoke a generic FDS simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "fdm.fds.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        post_processing_filename: str = None,
        n_cores: int = 1,
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[types.Path] = "",
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           input_filename=sim_config_filename,
                           post_processing_config=post_processing_filename,
                           storage_dir=storage_dir,
                           n_cores=n_cores)
