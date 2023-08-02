"""DualSPHysics module of the API."""
from typing import Literal, Optional
from uuid import UUID

from inductiva import types
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
        resource_pool_id: Optional[UUID] = None,
        device: Literal["gpu", "cpu"] = "cpu",
        run_async: bool = False,
    ):
        """Run the simulation.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            device: Device in which to run the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            resource_pool_id=resource_pool_id,
            device=device,
            input_filename=sim_config_filename,
            run_async=run_async,
        )
