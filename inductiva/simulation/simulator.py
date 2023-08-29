"""Base class for low-level simulators."""
from typing import Optional
from abc import ABC, abstractmethod

from inductiva import types, tasks, resources
from inductiva.utils import files


class Simulator(ABC):
    """Base simulator class."""

    @property
    @abstractmethod
    def api_method_name(self) -> str:
        pass

    def _setup_input_dir(self, input_dir: types.Path):
        """Setup the simulator input directory."""
        input_dir = files.resolve_path(input_dir)
        if not input_dir.is_dir():
            raise ValueError(
                f"The provided path (\"{input_dir}\") is not a directory.")
        return input_dir

    def run(
        self,
        input_dir: types.Path,
        *_args,
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        **kwargs,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory containing the input files.
            _args: Unused in this method, but defined to allow for more
                non-default arguments in method override in subclasses.
            **kwargs: Additional keyword arguments to be passed to the
                simulation API method.
        """
        input_dir = self._setup_input_dir(input_dir)

        return tasks.run_simulation(
            self.api_method_name,
            input_dir,
            run_async=run_async,
            machine_group=machine_group,
            **kwargs,
        )
