"""Base class for low-level simulators."""
from abc import ABC, abstractmethod
import pathlib
from typing import Optional

from inductiva import api
from inductiva import types
from inductiva.utils import files


class Simulator(ABC):
    """Base simulator class."""

    @property
    @abstractmethod
    def api_method_name(self) -> str:
        pass

    def run(
        self,
        input_dir: types.Path,
        *_args,  # unused in this method, but defined to allow for more
        # non-default arguments in method override in subclasses.
        output_dir: Optional[types.Path] = None,
        **kwargs,
    ) -> pathlib.Path:
        """Run the simulation."""
        input_dir = files.resolve_path(input_dir)
        if not input_dir.is_dir():
            raise ValueError(
                f"The provided path (\"{input_dir}\") is not a directory.")

        if output_dir is None:
            output_dir = input_dir.with_name(f"{input_dir.name}-output")
            output_dir = files.get_timestamped_path(output_dir)
        else:
            output_dir = files.resolve_path(output_dir)

        return api.run_simulation(
            self.api_method_name,
            input_dir,
            output_dir,
            params=kwargs,
        )
