"""SplisHSPlasH simulator module of the API."""
from typing import Literal, Optional

from inductiva import types, tasks, resources
from inductiva.simulators import Simulator


class SplishSplash(Simulator):
    """Class to invoke a generic SPlisHSPlasH simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sph.splishsplash.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        machine_group: Optional[resources.MachineGroup] = None,
        particle_radius: float = 0.025,
        device: Literal["gpu", "cpu"] = "cpu",
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            device: Device in which to run the simulation.
            particle_radius: Radius of the particles used in the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            machine_group=machine_group,
            particle_radius=particle_radius,
            device=device,
            input_filename=sim_config_filename,
            run_async=run_async,
        )
