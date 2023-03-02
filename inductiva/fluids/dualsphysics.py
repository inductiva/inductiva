"""DualSPHysics module of the API."""
import os
import pathlib

import inductiva
from inductiva.types import Path


class DualSPHysics:
    """Class to invoke a generic DualSPHysics simulation on the API.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        input_filename: Name of the DualSPHysics input file. The file should
            be present in `sim_dir`, and the name is relative to that
            directory. By default it is `DualSPHysics_input.xml`.
    """

    def __init__(
        self,
        sim_dir: Path,
        input_filename = "InputCase_Def.xml",
    ):
        self.input_filename = input_filename
        self.sim_dir = pathlib.Path(sim_dir)

        if not os.path.isdir(sim_dir):
            raise ValueError("The provided path is not a directory.")

    def simulate(self, output_dir=None) -> Path:
        """Run the simulation.

        Args:
            output_dir: Directory where the generated files will be stored.
        """
        return inductiva.sph.dualsphysics.run_simulation(self.sim_dir,
                                                         self.input_filename,
                                                         output_dir=output_dir)