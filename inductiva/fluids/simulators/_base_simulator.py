"""Base class for low-level simulators."""
import abc
import pathlib

from inductiva.types import Path


class BaseSimulator(abc.ABC):
    """Base class for the low-level simulator classes.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        input_filename: Name of the simulatoin input file. The file should
            be present in `sim_dir`, and the name is relative to that
            directory.
    """

    def __init__(self, sim_dir: Path, input_filename: str):
        self.sim_dir = pathlib.Path(sim_dir)
        self.input_filename = input_filename

        if not self.sim_dir.is_dir():
            raise ValueError(
                f"The provided path ({self.sim_dir}) is not a directory.")

        input_file_path = self.sim_dir.joinpath(input_filename)

        if not input_file_path.is_file():
            raise ValueError(f"\"{input_filename}\" input file not found "
                             f" in \"{self.sim_dir}\" simulation directory.")
