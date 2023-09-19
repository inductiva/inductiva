"""Heat sink scenario."""
from functools import singledispatchmethod
import os
import shutil
from typing import Optional

from inductiva import tasks, resources, simulators, scenarios, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "heat_sink")
OPENFOAM_TEMPLATE_SUBDIR = "openfoam"
FILES_SUBDIR = "files"
OPENFOAM_TEMPLATE_PARAMS_FILE_NAME = "parameters.jinja"
OPENFOAM_PARAMS_FILE_NAME = "parameters"
COMMANDS_FILE_NAME = "commands.json"


class HeatSink(scenarios.Scenario):
    """Heat sink scenario.

    This is a simulation scenario for a heat sink. A heat source is placed
    in a box where there is an air flow. A heat sink, placed on top of the
    source, is used to dissipate the heat via convection with the air flow.

    The heat source is modeled as a heater with a given power. The heat sink
    is a block of aluminum with thin fins on top. The air flow is introduced
    in the simulation via an inlet, where the air is at a fixed temperature.

    The scenario is simulated in a 3D box with dimensions 8 x 16 x 52 cm in x, y
    and z directions, respectively. The heat sink is a 4 x 3 x 6 cm block
    centered in (x, z) in the simulation box. The sink has 1 mm wide fins
    elongated along the z direction and separated by 2 mm centered in the x
    direction. The heat source is a 1 x 1 x 1 cm cube that sits under the heat
    sink. The air flow is injected in the simulation from the lower z boundary,
    and flows in the positive z direction.

    Schematic representations of the simulation scenario:

    - as seen from the side (zy plane): z points right, y points up, x points
      into the screen.
      _________________________________
      |                               |
      |                               |
      |                               |
      |  air flow ->                  |
      |           _________           |
      |           |       | heat      |
      |           |       | sink      |
      |           |       |           |
      |___________|_______|___________|
                    |___|
                 heat source

    - as seen from the the inlet of the air flow (xy plane): x points right, y
      points up, z points out of the screen.

      _________________________________
      |                               |
      |                               |
      |                               |
      |                               |
      |                               |
      |      fins | | | | | heat      |
      |           | | | | | sink      |
      |           |_|_|_|_|           |
      |           |       |           |
      |___________|_______|___________|
                    |___|
                 heat source

    The scenario can be simulated with OpenFOAM.
    """

    valid_simulators = [simulators.OpenFOAM]

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
        simulator: simulators.Simulator = simulators.OpenFOAM(),
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        simulation_time=300,
        output_time_step=10,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            simulation_time: The simulation time, in seconds.
            output_time_step: The time step to save the simulation results, in
              seconds.
            machine_group: The machine group to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        simulator.override_api_method_prefix("heat_sink")

        self.simulation_time = simulation_time
        self.output_time_step = output_time_step

        commands = self.get_commands()

        # TODO: Address the mpirun runs of HeatSink
        task = super().simulate(simulator,
                                machine_group=machine_group,
                                run_async=run_async,
                                n_cores=1,
                                commands=commands)

        return task

    def get_commands(self):
        """Returns the OpenFOAM commands for the simulation.
        """
        commands_file_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                          OPENFOAM_TEMPLATE_SUBDIR,
                                          COMMANDS_FILE_NAME)

        commands = self.read_commands_from_file(commands_file_path)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: simulators.Simulator):
        pass


@HeatSink.create_input_files.register
def _(self, simulator: simulators.OpenFOAM, input_dir):  # pylint: disable=unused-argument
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

    utils.templates.replace_params(
        template_path=template_file_path,
        params={
            "simulation_time": self.simulation_time,
            "output_time_step": self.output_time_step,
            "air_velocity": self.air_velocity,
            "air_temperature": self.air_temperature,
            "heater_power": self.heater_power
        },
        output_file=params_file_path,
        remove_template=True,
    )
