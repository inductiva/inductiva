"""Heat sink scenario."""
from functools import singledispatchmethod
import os
import shutil
import tempfile
from typing import Optional

from inductiva.types import Path
from inductiva.fluids.simulators import OpenFOAM
from inductiva.simulation import Simulator
from inductiva.scenarios import Scenario
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       replace_params_in_template)

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "heat_sink")
OPENFOAM_TEMPLATE_SUBDIR = "openfoam"
FILES_SUBDIR = "files"
COMMANDS_TEMPLATE_FILE = "commands.json.jinja"


class HeatSink(Scenario):
    """Heat sink scenario.
    
    This is a simulation scenario for a heat sink. A heat source is placed
    in a box where there is an air flow. A heat sink, placed on top of the 
    source, is used to dissipate the heat via convection with the air flow.

    The heat source is modeled as a heater with a given power. The heat sink
    is a block of aluminum with thin fins on top. The air flow is introduced
    in the simulation via an inlet, where the air is at a fixed temperature.
    """

    def __init__(
        self,
        temperature=290,
        heater_power=200,
    ):
        """Initializes the heat sink scenario.
        Args:
            temperature: The temperature of the air flow, in Kelvin. Also sets
              the initial temperature of the heater and the heat sink.
            heater_power: The power of the heater, in Watts.
        """
        self.temperature = temperature
        self.heater_power = heater_power

    def simulate(
            self,
            simulator: Simulator = OpenFOAM(),
            output_dir: Optional[Path] = None,
    ):
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            output_dir: The output directory to save the simulation results.
        """
        commands = self.get_commands()

        return super().simulate(simulator, output_dir, commands=commands)

    def simulate_async(
            self,
            simulator: Simulator = OpenFOAM(),
    ):
        """Simulates the scenario asynchronously.

        Args:
            simulator: The simulator to use for the simulation.
            output_dir: The output directory to save the simulation results.
        """

        commands = self.get_commands()

        return super().simulate_async(simulator, commands=commands)

    def get_commands(self):
        """Returns the commands for the simulation.
        """

        openfoam_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              OPENFOAM_TEMPLATE_SUBDIR)

        # Replace parameters in the template file, writing the result to a
        # temporary file.
        with tempfile.NamedTemporaryFile(mode="w") as commands_file:
            replace_params_in_template(openfoam_template_path,
                                       COMMANDS_TEMPLATE_FILE, {
                                           "temperature": self.temperature,
                                           "heater_power": self.heater_power,
                                       }, commands_file.name)

            # Read the commands from the temporary file.
            commands = self.read_commands_from_file(commands_file.name)

        return commands

    @singledispatchmethod
    def gen_config(self, simulator: Simulator):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def get_config_filename(self, simulator: Simulator):  # pylint: disable=unused-argument
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )


@HeatSink.get_config_filename.register
def _(self, simulator: OpenFOAM):  # pylint: disable=unused-argument
    pass


@HeatSink.gen_aux_files.register
def _(self, simulator: OpenFOAM, input_dir):  # pylint: disable=unused-argument
    """Setup the working directory for the simulation."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_SUBDIR, FILES_SUBDIR)

    shutil.copytree(template_files_dir, input_dir, dirs_exist_ok=True)


@HeatSink.gen_config.register
def _(self, simulator: OpenFOAM, input_dir):  # pylint: disable=unused-argument
    """Generate the mdp configuration files for the simulation."""
    pass
