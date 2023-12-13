"""Base class for scenarios."""

from abc import ABC
import io
import json
import shutil
import os
import random
import string
from typing import Optional

import inductiva


class Scenario(ABC):
    """Base class for scenarios."""
    params = {}
    template_files_dir = None
    input_base_dir = "scenario-input-dir"

    def _create_input_dir(self):
        """Create input directory."""

    def pre_simulate_hook(self):
        """Hook to be executed before simulation.
        
        This method encapsulates all steps requried to be executed
        before the simulation is run.
        """
        pass

    def add_extra_input_files(self, input_paths=None, output_paths=None):
        """Entry-point to add extra files used in the simulation."""

        # Wrap-up in lists
        if not isinstance(input_paths, list):
            input_paths = [input_paths]
        if not isinstance(output_paths, list):
            output_paths = [output_paths]

        output_paths = [
            os.path.join(self.input_dir, path) for path in output_paths
        ]
        inductiva.utils.files.copy_files(input_paths, output_paths)

    def _create_input_files(self):
        """Create input files from template."""

        # Copy all files from the template dir to the input directory
        shutil.copytree(self.template_files_dir,
                        self.input_dir,
                        dirs_exist_ok=True,
                        symlinks=True)

        template_filenames = inductiva.utils.templates.get_template_filenames(
            self.template_files_dir)

        output_filename_paths = [
            os.path.join(self.input_dir,
                         file.split(".jinja")[0]) for file in template_filenames
        ]

        inductiva.utils.templates.batch_replace_params(
            templates_dir=self.input_dir,
            template_filenames=template_filenames,
            params=self.params,
            output_filename_paths=output_filename_paths,
            remove_templates=True,
        )

    def _create_command_file(self):
        """Create command file from template."""

        commands_file = os.path.join(self.template_files_dir, "commands.json")
        commands_template_path = commands_file + ".jinja"

        if os.path.isfile(commands_template_path):
            inmemory_file = io.StringIO()
            inductiva.utils.templates.replace_params(
                template_path=commands_template_path,
                params=self.params,
                output_file=inmemory_file,
            )
            return inmemory_file

        return commands_file

    def get_commands(self):
        "Read list of commands from commands.json file"

        commands_file = self._create_command_file()

        if isinstance(commands_file, str):
            if os.path.exists(commands_file):
                with open(commands_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            else:
                return None
        # Make sure already opened file is read from the beginning
        commands_file.seek(0)
        return json.load(commands_file)

    def simulate(
        self,
        machine_group: Optional[inductiva.resources.MachineGroup] = None,
        storage_dir: Optional[inductiva.types.Path] = "",
        **kwargs,
    ):
        """Run simulation."""

        # Create input directory
        self.input_dir = \
            f"{self.input_base_dir}-{inductiva.utils.misc.create_random_tag()}"
        os.mkdir(self.input_dir)

        self.pre_simulate_hook()
        self._create_input_files()

        return self.simulator.run(
            input_dir=self.input_dir,
            machine_group=machine_group,
            storage_dir=storage_dir,
            **kwargs,
        )
