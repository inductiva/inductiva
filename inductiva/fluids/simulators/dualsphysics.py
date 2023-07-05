"""DualSPHysics module of the API."""
import pathlib
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
        output_dir: Optional[types.Path] = None,
        device: Literal["gpu", "cpu"] = "cpu",
        resource_pool_id: Optional[UUID] = None,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            device: Device in which to run the simulation.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           output_dir=output_dir,
                           resource_pool_id=resource_pool_id,
                           device=device,
                           input_filename=sim_config_filename)

    def run_async(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        device: Literal["gpu", "cpu"] = "cpu",
        resource_pool_id: Optional[UUID] = None,
    ) -> str:
        """Run the simulation asynchronously.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            device: Device in which to run the simulation.
            other arguments: See the documentation of the base class.
            """

        return super().run_async(input_dir,
                                 resource_pool_id=resource_pool_id,
                                 device=device,
                                 input_filename=sim_config_filename)
