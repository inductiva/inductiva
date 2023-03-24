"""DualSPHysics module of the API."""
import pathlib
from typing import Literal, Optional

from inductiva.types import Path
from inductiva.fluids.simulators._base_simulator import BaseSimulator


class SPlisHSPlasH(BaseSimulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self, sim_dir: Path, sim_config_filename: str):
        super().__init__(sim_dir, sim_config_filename,
                         "sph.splishsplash.run_simulation")

    def simulate(
        self,
        output_dir: Optional[Path] = None,
        device: Literal["gpu", "cpu"] = "cpu",
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            device: Device in which to run the simulation."""
        return super().simulate(output_dir=output_dir, device=device)
