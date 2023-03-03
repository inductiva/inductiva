"""SWASH module of the API."""
import os
import pathlib

import inductiva
from inductiva.types import Path


class SWASH:
    """Class to invoke a generic SWASH simulation on the API.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        input_filename: Name of the SWASH input file. The file should
            be present in `sim_dir`, and the name is relative to that
            directory.
    """

    def __init__(
        self,
        sim_dir: Path,
        input_filename: str,
    ):
        self.input_filename = input_filename
        self.sim_dir = pathlib.Path(sim_dir)

        if not os.path.isdir(sim_dir):
            raise ValueError("The provided path is not a directory.")

    def simulate(self, n_cores=1, output_dir=None) -> Path:
        """Run the simulation.

        Args:
            output_dir: Directory where the generated files will be stored.
        """
        return inductiva.shallow.swash.run_simulation(self.sim_dir,
                                                      self.input_filename,
                                                      n_cores=n_cores,
                                                      output_dir=output_dir)
