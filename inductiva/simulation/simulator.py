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
        *_args,
        output_dir: Optional[types.Path] = None,
        track_logs: bool = False,
        **kwargs,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            input_dir: Path to the directory containing the input files.
            _args: Unused in this method, but defined to allow for more
                non-default arguments in method override in subclasses.
            output_dir: Path to the directory where the output files will be
                stored. If not provided, a timestamped directory will be
                created with the same name as the input directory appended
                with "-output".
            track_logs: If True, the logs of the remote execution will be
                streamed to the console.
            **kwargs: Additional keyword arguments to be passed to the
                simulation API method.
        """
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
            log_remote_execution=track_logs,
        )

    def run_async(
        self,
        input_dir: types.Path,
        *_args,
        **kwargs,
    ) -> str:
        """Run the simulation asynchronously.

        Args:
            input_dir: Path to the directory containing the input files.
            _args: Unused in this method, but defined to allow for more
                non-default arguments in method override in subclasses.
            **kwargs: Additional keyword arguments to be passed to the
                simulation API method.
        """
        input_dir = files.resolve_path(input_dir)
        if not input_dir.is_dir():
            raise ValueError(
                f"The provided path (\"{input_dir}\") is not a directory.")

        return api.run_async_simulation(
            self.api_method_name,
            input_dir,
            params=kwargs,
        )
