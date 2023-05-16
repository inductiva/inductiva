"""DualSPHysics module of the API."""
import pathlib
from typing import Literal, Optional

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
        output_dir: Optional[types.Path] = None,
        device: Literal["gpu", "cpu"] = "cpu",
        track_logs: bool = False,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            device: Device in which to run the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            output_dir=output_dir,
            track_logs=track_logs,
            device=device,
            input_filename=sim_config_filename,
        )

    def run_async(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        device: Literal["gpu", "cpu"] = "cpu",
    ) -> str:
        """Run the simulation asynchronously.
        
        Args:
            sim_config_filename: Name of the simulation configuration file.
            device: Device in which to run the simulation.
            other arguments: See the documentation of the base class.
            """

        return super().run_async(input_dir,
                                 device=device,
                                 input_filename=sim_config_filename)
