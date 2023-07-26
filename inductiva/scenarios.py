"""Base class for scenarios."""

from abc import ABC, abstractmethod
import tempfile
from typing import Optional, Union
from uuid import UUID
from inductiva.types import Path
from inductiva.simulation import Simulator
from inductiva.utils.misc import split_camel_case
import json
from inductiva.tasks import Task


class Scenario(ABC):
    """Base class for scenarios."""
    valid_simulators = []

    def _output_dir_name(self, output_dir: str):
        """Setup the scenario output directory."""
        if output_dir is None:
            scenario_name_splitted = split_camel_case(self.__class__.__name__)
            return "-".join(scenario_name_splitted).lower()
        else:
            return output_dir

    @abstractmethod
    def create_input_files(self, simulator: Simulator, input_dir: Path):
        """To be implemented in subclasses."""
        pass

    def read_commands_from_file(self, commands_path: str):
        "Read list of commands from commands.json file"
        with open(commands_path, "r", encoding="utf-8") as f:
            commands = json.load(f)
        return commands

    def validate_simulator(self, simulator: Simulator):
        """Checks if the scenario can be simulated with the given simulator."""
        if type(simulator) not in self.valid_simulators:
            raise ValueError(
                f"Simulator not supported for `{self.__class__.__name__}` "
                "scenario.")

    def simulate(
        self,
        simulator: Simulator,
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
        **kwargs,
    ) -> Union[Path, Task]:
        """Simulates the scenario synchronously."""
        self.validate_simulator(simulator)

        with tempfile.TemporaryDirectory() as input_dir:
            if not run_async:
                self.create_input_files(simulator, input_dir)
                output_dir = self._output_dir_name(output_dir)

            return simulator.run(
                input_dir,
                output_dir=output_dir,
                resource_pool_id=resource_pool_id,
                run_async=run_async,
                **kwargs,
            )
