"""Heat sink scenario."""
from functools import singledispatchmethod
import os
import shutil
from typing import Optional
from uuid import UUID

from inductiva.tasks import Task
from inductiva.types import Path
from inductiva.fluids.simulators import OpenFOAM
from inductiva.simulation import Simulator
from inductiva.scenarios import Scenario
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       replace_params_in_template)

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "heat_sink")
OPENFOAM_TEMPLATE_SUBDIR = "openfoam"
FILES_SUBDIR = "files"
OPENFOAM_TEMPLATE_PARAMS_FILE_NAME = "parameters.jinja"
OPENFOAM_PARAMS_FILE_NAME = "parameters"
COMMANDS_FILE_NAME = "commands.json"


class HeatSink(Scenario):
    """Heat sink scenario.

    This is a simulation scenario for a heat sink. A heat source is placed
    in a box where there is an air flow. A heat sink, placed on top of the
    source, is used to dissipate the heat via convection with the air flow.

    The heat source is modeled as a heater with a given power. The heat sink
    is a block of aluminum with thin fins on top. The air flow is introduced
    in the simulation via an inlet, where the air is at a fixed temperature.
    """

    valid_simulators = [OpenFOAM]

    def __init__(
        self,
        air_velocity=10,
        air_temperature=290,
        heater_power=200,
    ):
        """Initializes the heat sink scenario.

        Args:
            air_velocity: The velocity of the air flow, in m/s.
            air_temperature: The temperature of the air flow, in Kelvin. Also
              sets the initial temperature of the heater and the heat sink.
            heater_power: The power of the heater, in Watts.
        """
        self.air_velocity = air_velocity
        self.air_temperature = air_temperature
        self.heater_power = heater_power

    def simulate(
        self,
        simulator: Simulator = OpenFOAM(),
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        simulation_time=300,
        output_time_step=10,
    ):
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            output_dir: The output directory to save the simulation results.
            simulation_time: The simulation time, in seconds.
            output_time_step: The time step to save the simulation results, in
              seconds.
        """
        self.simulation_time = simulation_time
        self.output_time_step = output_time_step

        commands = self.get_commands()

        return super().simulate(simulator,
                                output_dir,
                                resource_pool_id=resource_pool_id,
                                commands=commands)

    def simulate_async(
        self,
        simulator: Simulator = OpenFOAM(),
        resource_pool_id: Optional[UUID] = None,
        simulation_time=300,
        output_time_step=10,
    ) -> Task:
        """Simulates the scenario asynchronously.

        Args:
            simulator: The simulator to use for the simulation.
            output_dir: The output directory to save the simulation results.
            simulation_time: The simulation time, in seconds.
            output_time_step: The time step to save the simulation results, in
              seconds.
        """
        self.simulation_time = simulation_time
        self.output_time_step = output_time_step

        commands = self.get_commands()

        return super().simulate_async(simulator,
                                      resource_pool_id=resource_pool_id,
                                      commands=commands)

    def get_commands(self):
        """Returns the OpenFOAM commands for the simulation.
        """
        commands_file_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                          OPENFOAM_TEMPLATE_SUBDIR,
                                          COMMANDS_FILE_NAME)

        commands = self.read_commands_from_file(commands_file_path)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@HeatSink.create_input_files.register
def _(self, simulator: OpenFOAM, input_dir):  # pylint: disable=unused-argument
    """Creates OpenFOAM simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_SUBDIR, FILES_SUBDIR)

    shutil.copytree(template_files_dir,
                    input_dir,
                    dirs_exist_ok=True,
                    symlinks=True)

    template_file_path, params_file_path = (
        os.path.join(input_dir, file_name) for file_name in
        [OPENFOAM_TEMPLATE_PARAMS_FILE_NAME, OPENFOAM_PARAMS_FILE_NAME])

    replace_params_in_template(
        template_path=template_file_path,
        params={
            "simulation_time": self.simulation_time,
            "output_time_step": self.output_time_step,
            "air_velocity": self.air_velocity,
            "air_temperature": self.air_temperature,
            "heater_power": self.heater_power
        },
        output_file_path=params_file_path,
        remove_template=True,
    )
