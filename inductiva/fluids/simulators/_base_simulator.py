"""Base class for low-level simulators."""
import abc
import pathlib
from typing import Optional

from inductiva_web_api_client.models import TaskRequest

from inductiva.api.methods import invoke_api
from inductiva.types import Path
from inductiva.utils.data import pack_input


class BaseSimulator(abc.ABC):
    """Base class for the low-level simulator classes.

    Attributes:
        sim_dir: Path to the directory with all the simulation input files.
        main_config_filename: Name of the main simulation config file. The file
            should be present in `sim_dir`, and the name is relative to that
            directory.
    """

    @abc.abstractmethod
    def __init__(
        self,
        sim_dir: Path,
        sim_config_filename: str,
        api_method_name: str,
    ):
        self.sim_dir = pathlib.Path(sim_dir)
        self.sim_config_filename = sim_config_filename
        self.api_method_name = api_method_name

        if not self.sim_dir.is_dir():
            raise ValueError(
                f"The provided path ({self.sim_dir}) is not a directory.")

        input_file_path = self.sim_dir.joinpath(sim_config_filename)

        if not input_file_path.is_file():
            raise ValueError(f"\"{sim_config_filename}\" input file not found "
                             f" in \"{self.sim_dir}\" simulation directory.")

    def simulate(self, output_dir: Optional[Path] = None, **kwargs):
        params = {
            "sim_dir": str(self.sim_dir),
            "input_filename": self.sim_config_filename,
            **kwargs
        }

        request = TaskRequest(method=self.api_method_name, params=params)

        input_zip = pack_input(params, {"sim_dir": Path})

        return invoke_api(request, input_zip, output_dir=output_dir)
