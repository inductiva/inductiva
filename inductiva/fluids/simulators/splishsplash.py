"""DualSPHysics module of the API."""
from typing import Literal, Optional

from inductiva import types, tasks, resources
from inductiva.simulation import Simulator


class SPlisHSPlasH(Simulator):
    """Class to invoke a generic SPlisHSPlasH simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sph.splishsplash.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        machine_group: Optional[resources.MachineGroup] = None,
        device: Literal["gpu", "cpu"] = "cpu",
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            device: Device in which to run the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            machine_group=machine_group,
            device=device,
            input_filename=sim_config_filename,
            run_async=run_async,
        )
