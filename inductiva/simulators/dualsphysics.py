"""DualSPHysics simulator module of the API."""

from typing import Optional

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
        sim_config_filename: str = "config.xml",
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
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
                           run_async=run_async)
