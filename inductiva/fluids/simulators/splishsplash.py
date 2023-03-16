"""SplishSplash module of the API."""
import os
import pathlib
from typing import Literal

import inductiva
from inductiva.types import Path


class SPlisHSPlasH:
    """Class to invoke a generic SPlisHSPlasH simulation on the API.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        input_filename: Name of the SPlisHSPlasH input file. The file should
            be present in `sim_dir`, and the name is relative to that
            directory. By default it is `splishsplash_input.json`.
    """

    def __init__(
        self,
        sim_dir: Path,
        input_filename: str = "splishsplash_input.json",
    ):
        self.input_filename = input_filename
        self.sim_dir = pathlib.Path(sim_dir)

        if not os.path.isdir(sim_dir):
            raise ValueError("The provided path is not a directory.")

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
