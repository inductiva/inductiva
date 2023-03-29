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
        sim_dir: types.Path,
        output_dir: Optional[types.Path] = None,
        device: Literal["gpu", "cpu"] = "cpu",
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            device: Device in which to run the simulation.
        """
        return super().run(sim_dir, output_dir=output_dir, device=device)
