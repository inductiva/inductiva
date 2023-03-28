"""Base class for low-level simulators."""
from abc import ABC, abstractmethod
import pathlib
from typing import Optional

from inductiva.api.methods import invoke_api
from inductiva import types
from inductiva.utils.files import get_timestamped_path


class Simulator(ABC):
    """Base class for the low-level simulator classes.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        main_config_filename: Name of the main simulation config file. The file
            should be present in `sim_dir`, and the name is relative to that
            directory.
        api_method_name: Name of the API method that will be invoked to run the
            simulation.
    """

    def __init__(
        self,
        sim_dir: types.Path,
        sim_config_filename: str,
    ):
        self.sim_dir = pathlib.Path(sim_dir)
        self.sim_config_filename = sim_config_filename

        if not self.sim_dir.is_dir():
            raise ValueError(
                f"The provided path ({self.sim_dir}) is not a directory.")

        input_file_path = self.sim_dir.joinpath(sim_config_filename)

        if not input_file_path.is_file():
            raise ValueError(f"\"{sim_config_filename}\" input file not found "
                             f" in \"{self.sim_dir}\" simulation directory.")

    @property
    @abstractmethod
    def api_method_name(self) -> str:
        pass

    def simulate(self,
                 output_dir: Optional[types.Path] = None,
                 **kwargs) -> pathlib.Path:
        params = {
            "sim_dir": self.sim_dir,
            "input_filename": self.sim_config_filename,
            **kwargs
        }

        type_annotations = {"sim_dir": types.Path}

        if output_dir is None:
            sim_dir_name = self.sim_dir.name
            output_dir = self.sim_dir.with_name(f"{sim_dir_name}-output")
            output_dir = get_timestamped_path(output_dir)

        return invoke_api(self.api_method_name,
                          params,
                          type_annotations,
                          output_dir=output_dir)
