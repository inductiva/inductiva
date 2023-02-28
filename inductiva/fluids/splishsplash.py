"""SplishSplash module of the API."""
import os
import pathlib

import inductiva
from inductiva.types import Path


class SPlisHSPlasH:
    """Class to invoke a generic SplishSplash simulation on the API."""

    def __init__(self, sim_dir: Path, input_file: str):
        self.input_file = input_file
        self.sim_dir = pathlib.Path(sim_dir)

        if not os.path.isdir(sim_dir):
            raise ValueError("The provided path is not a directory.")

    def simulate(self, output_dir=None) -> Path:
        return inductiva.sph.splishsplash.run_simulation(self.sim_dir,
                                                         self.input_file,
                                                         output_dir=output_dir)
