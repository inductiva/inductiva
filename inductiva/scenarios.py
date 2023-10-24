"""Base class for scenarios."""

from abc import ABC, abstractmethod
import io
import shutil
import glob
import os
import tempfile
from typing import Optional, Union

from inductiva import resources, utils
from inductiva.types import Path
from inductiva.simulators import Simulator
import json


class Scenario(ABC):
    """Base class for scenarios."""
    valid_simulators = []
    params = {}
    template_files_dir = None

    def config_input(self, simulator: Simulator, input_dir: Path):
        """Possible to be implemented in subclasses."""
        pass

    def add_input_files(self, simulator: Simulator, input_dir: Path):
        """Possible to be implemented in subclasses."""
        pass

    def create_input_files(self, simulator: Simulator, input_dir: Path):
        """Create input files from template."""

        template_files_dir = os.path.join(self.template_files_dir, "sim_config_files")

        # Copy all files from the template dir to the input directory
        shutil.copytree(template_files_dir,
                        input_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        template_filenames = glob.glob(
            os.path.join("**", "*.jinja"),
            root_dir=template_files_dir,
            recursive=True)
        output_filename_paths = [
            os.path.join(input_dir, file.split(".jinja")[0]) for file in template_filenames
        ]
        utils.templates.batch_replace_params(
            templates_dir=input_dir,
            template_filenames=template_filenames,
            params=self.params,
            output_filename_paths=output_filename_paths,
            remove_templates=True,
        )

    def get_commands(self, commands_file: Union[str, io.StringIO] = "commands.json"):
        "Read list of commands from commands.json file"

        commands_file = os.path.join(self.template_files_dir, commands_file)
        if isinstance(commands_file, str):
            with open(commands_file, "r", encoding="utf-8") as f:
                return json.load(f)

        # Make sure already opened file is read from the beginning
        commands_file.seek(0)
        return json.load(commands_file)

    def validate_simulator(self, simulator: Simulator):
        """Checks if the scenario can be simulated with the given simulator."""
        if type(simulator) not in self.valid_simulators:
            raise ValueError(
                f"Simulator not supported for `{self.__class__.__name__}` "
                "scenario.")

    def simulate(
        self,
        simulator: Simulator,
        machine_group: Optional[resources.MachineGroup] = None,
        **kwargs,
    ):
        """Simulates the scenario synchronously."""
        self.validate_simulator(simulator)

        with tempfile.TemporaryDirectory() as input_dir:
            self.config_input(simulator, input_dir)
            self.create_input_files(simulator, input_dir)
            self.add_input_files(simulator, input_dir)

            return simulator.run(
                input_dir,
                machine_group=machine_group,
                **kwargs,
            )
