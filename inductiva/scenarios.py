"""Base class for scenarios."""

from abc import ABC
import io
import json
import shutil
import os
import tempfile
from typing import Optional

from inductiva import simulators, resources, types, utils


class Scenario(ABC):
    """Base class for scenarios."""
    valid_simulators = []
    params = {}
    template_files_dir = None

    def config_params(self, simulator: simulators.Simulator,
                      input_dir: types.Path):
        """Entry-point to further configure params.
        
        Useful when a scenario can be configured for several simulators,
        but params differ in implementation between them.
        """
        pass

    def add_extra_input_files(self, simulator: simulators.Simulator,
                              input_dir: types.Path):
        """Entry-point to add extra files used in the simulation.
        
        Usefull to files as args to the simulation. E.g., protein or vehicle.
        """
        pass

    def create_input_files(
            self,
            simulator: simulators.Simulator,  # pylint: disable=unused-argument
            input_dir: types.Path):
        """Create input files from template."""

        template_files_dir = os.path.join(self.template_files_dir,
                                          "sim_config_files")

        # Copy all files from the template dir to the input directory
        shutil.copytree(template_files_dir,
                        input_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        template_filenames = utils.templates.get_template_filenames(
            template_files_dir)

        output_filename_paths = [
            os.path.join(input_dir,
                         file.split(".jinja")[0]) for file in template_filenames
        ]

        utils.templates.batch_replace_params(
            templates_dir=input_dir,
            template_filenames=template_filenames,
            params=self.params,
            output_filename_paths=output_filename_paths,
            remove_templates=True,
        )

    def create_command_file(self):
        """Create command file from template."""

        commands_file = os.path.join(self.template_files_dir, "commands.json")
        commands_template_path = commands_file + ".jinja"

        if os.path.isfile(commands_template_path):
            inmemory_file = io.StringIO()
            utils.templates.replace_params(
                template_path=commands_template_path,
                params=self.params,
                output_file=inmemory_file,
            )
            return inmemory_file

        return commands_file

    def get_commands(self):
        "Read list of commands from commands.json file"

        commands_file = self.create_command_file()

        if isinstance(commands_file, str):
            if os.path.exists(commands_file):
                with open(commands_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            else:
                return None
        # Make sure already opened file is read from the beginning
        commands_file.seek(0)
        return json.load(commands_file)

    def validate_simulator(self, simulator: simulators.Simulator):
        """Checks if the scenario can be simulated with the given simulator."""
        if type(simulator) not in self.valid_simulators:
            raise ValueError(
                f"Simulator not supported for `{self.__class__.__name__}` "
                "scenario.")

    def simulate(
        self,
        simulator: simulators.Simulator,
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[types.Path] = "",
        **kwargs,
    ):
        """Simulates the scenario synchronously."""
        self.validate_simulator(simulator)

        with tempfile.TemporaryDirectory() as input_dir:
            self.config_params(simulator, input_dir)
            self.create_input_files(simulator, input_dir)
            self.add_extra_input_files(simulator, input_dir)

            kwargs = {
                key: value for key, value in kwargs.items() if value is not None
            }

            return simulator.run(
                input_dir,
                machine_group=machine_group,
                storage_dir=storage_dir,
                **kwargs,
            )
