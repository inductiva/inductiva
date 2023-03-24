"""SplishSplash module of the API."""
from typing import Literal

import inductiva
from inductiva.fluids.simulators._base_simulator import BaseSimulator
from inductiva.types import Path


class SPlisHSPlasH(BaseSimulator):
    """Class to invoke a generic SPlisHSPlasH simulation on the API."""

    def __init__(
        self,
        sim_dir: Path,
        input_filename: str = "splishsplash_input.json",
    ):
        super().__init__(sim_dir, input_filename)

    def simulate(self,
                 device: Literal["gpu", "cpu"] = "cpu",
                 output_dir=None) -> Path:
        """Run the simulation.

        Args:
            output_dir: Directory where the generated files will be stored.
        """
        return inductiva.sph.splishsplash.run_simulation(self.sim_dir,
                                                         self.input_filename,
                                                         device=device,
                                                         output_dir=output_dir)
