"""Base class for scenarios."""

from abc import ABC, abstractmethod
import tempfile
from typing import Optional
from uuid import UUID
from inductiva.types import Path
from inductiva.simulation import Simulator
from inductiva.utils.files import resolve_path, get_timestamped_path
from inductiva.utils.misc import split_camel_case
from inductiva.tasks.methods import get_task_info, fetch_task_output
import json


class Scenario(ABC):
    """Base class for scenarios."""

    def _setup_output_dir(self, output_dir: str):
        """Setup the scenario output directory."""
        if output_dir is None:
            scenario_name_splitted = split_camel_case(self.__class__.__name__)
            output_dir_prefix = "-".join(scenario_name_splitted).lower()
            output_dir = get_timestamped_path(f"{output_dir_prefix}-output")
        output_dir = resolve_path(output_dir)
        return output_dir

    def _setup_config(self, simulator: Simulator, input_dir: Path):
        """Setup the scenario configuration files and arguments."""
        self.gen_aux_files(simulator, input_dir)
        self.gen_config(simulator, input_dir)

        args = ()
        config_filename = self.get_config_filename(simulator)
        if config_filename:
            args += (config_filename,)
        return args

    def download_outputs(self, output_dir: str = None):
        """Download the outputs of an async simulation to output_dir."""

        if not get_task_info(self.task_id)["status"]:
            raise ValueError("Simulation not finished.")

        fetch_task_output(self.task_id, output_dir=output_dir, return_type=None)

    @abstractmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: Path):
        """To be implemented in subclasses."""
        pass

    @abstractmethod
    def gen_config(self, simulator: Simulator, input_dir: Path):
        """To be implemented in subclasses."""
        pass

    @abstractmethod
    def get_config_filename(self, simulator: Simulator):
        """To be implemented in subclasses."""
        pass

    def read_commands_from_file(self, commands_path: str):
        "Read list of commands from commands.json file"
        with open(commands_path, "r", encoding="utf-8") as f:
            commands = json.load(f)
        return commands

    def simulate(
        self,
        simulator: Simulator,
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        **kwargs,
    ) -> Path:
        """Simulates the scenario for a single simulator call."""
        output_dir = self._setup_output_dir(output_dir)
        with tempfile.TemporaryDirectory() as input_dir:
            args = self._setup_config(simulator, input_dir)
            return simulator.run(
                input_dir,
                *args,
                output_dir=output_dir,
                resource_pool_id=resource_pool_id,
                **kwargs,
            )

    def simulate_async(
        self,
        simulator: Simulator,
        resource_pool_id: Optional[UUID] = None,
        **kwargs,
    ) -> str:
        """Simulates the scenario asychronously."""

        with tempfile.TemporaryDirectory() as input_dir:
            args = self._setup_config(simulator, input_dir)
            task_id = simulator.run_async(
                input_dir,
                *args,
                resource_pool_id=resource_pool_id,
                **kwargs,
            )
        return task_id
