"""DualSPHysics module of the API."""
from typing import Optional

from inductiva import types, tasks, resources
from inductiva.simulators import Simulator


class XBeach(Simulator):
    """Class to invoke a generic XBeach simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sw.xbeach.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: Optional[str] = "params.txt",
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            input_filename=sim_config_filename,
            machine_group=machine_group,
            run_async=run_async,
        )
