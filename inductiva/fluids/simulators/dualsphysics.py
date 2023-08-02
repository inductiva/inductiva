"""DualSPHysics module of the API."""

from typing import Literal, Optional
from uuid import UUID

from inductiva import types
from inductiva.simulation import Simulator


class DualSPHysics(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sph.dualsphysics.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        device: Literal["gpu", "cpu"] = "cpu",
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
    ):
        """Run the simulation.

        Args:
            device: Device in which to run the simulation.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           resource_pool_id=resource_pool_id,
                           device=device,
                           input_filename=sim_config_filename,
                           run_async=run_async)
