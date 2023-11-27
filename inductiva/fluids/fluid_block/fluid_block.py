"""Fluid Block scenario."""
from functools import singledispatchmethod
import os
from typing import List, Optional
import shutil

from inductiva import tasks, resources, fluids, scenarios, simulators, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "fluid_block")
SPLISHPLASH_TEMPLATE_INPUT_DIR = "splishsplash"
SPLISHSPLASH_CONFIG_FILENAME = "fluid_block.json"
DUALSPHYSICS_TEMPLATE_INPUT_DIR = "dualsphysics"
DUALSPHYSICS_CONFIG_FILENAME = "fluid_block.xml"


class FluidBlock(scenarios.Scenario):
    """Fluid block scenario.

    This is a simulation scenario for a fluid block moving in a cubic tank under
    the action of gravity. The tank is a cube of dimensions 1 x 1 x 1 m. The
    fluid block is also cubic, but has configurable dimensions and initial
    position and velocity. The fluid properties such as density and kinematic
    viscosity are also configurable.

    Schematic representation of the simulation scenario:
    _________________________________
    |                               |  tank
    |         ___________           |
    |         |         |           |
    |         |  fluid  |  ->       |
    |         |  block  |  initial  |
    |         |_________|  velocity |
    |                               |
    |                               |
    |_______________________________|

    The scenario can be simulated with SPlisHSPlasH and DualSPHysics.
    """

    valid_simulators = [simulators.SplishSplash, simulators.DualSPHysics]

    def __init__(self,
                 density: float,
                 kinematic_viscosity: float,
                 dimensions: List[float],
                 position: Optional[List[float]] = None,
                 initial_velocity: Optional[List[float]] = None):
        """Initializes the fluid block scenario.

        Args:
            density: The density of the fluid in kg/m^3. Valid
                range is [400, 2000] Kg/m3.
            kinematic_viscosity: The kinematic viscosity of the fluid, in m^2/s.
                Reference value for water is 1e-6 m^2/s.
            dimensions: The fluid block dimensions (in x, y, z), in meters.
            position: The position of the fluid block in the tank (in x, y, z),
              in meters.
            initial_velocity: The initial velocity of the fluid block (in x, y,
              z), in m/s.
        """

        if not 400 <= density <= 2000:
            raise ValueError("`density` must be in range [400, 2000] Kg/m3.")

        self.params["fluid"] = fluids.FluidType(
            density=density, kinematic_viscosity=kinematic_viscosity)

        if len(dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")
        self.params["fluid_dimensions"] = dimensions

        if position is None:
            position = [0, 0, 0]
        self.params["fluid_position"] = position

        if initial_velocity is None:
            initial_velocity = [0, 0, 0]
        self.params["initial_velocity"] = initial_velocity

    def simulate(
        self,
        simulator: simulators.Simulator = simulators.DualSPHysics(),
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[str] = "",
        particle_radius: float = 0.02,
        simulation_time: float = 1,
        adaptive_time_step: bool = True,
        particle_sorting: bool = True,
        time_step: float = 0.001,
        output_export_rate: float = 60,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation. Supported
              simulators are: SPlisHSPlasH, DualSPHysics.
            machine_group: The machine group to use for the simulation.
            particle_radius: Radius of the fluid particles, in meters.
              Determines the resolution of the simulation. Lower values result
              in higher resolution and longer simulation times.
            simulation_time: The simulation time, in seconds.
            adaptive_time_step: Whether to use adaptive time stepping.
            particle_sorting: Whether to use particle sorting.
            time_step: Time step, in seconds.
            output_export_rate: Rate to export outputs per second.
            storage_dir: The parent directory where the simulation
            results will be stored.

        """
        simulator.override_api_method_prefix("fluid_block")

        self.params["particle_radius"] = particle_radius
        self.params["simulation_time"] = simulation_time
        self.params["adaptive_time_step"] = adaptive_time_step
        self.params["particle_sorting"] = particle_sorting
        self.params["time_step"] = time_step
        self.params["output_export_rate"] = output_export_rate

        self.set_template_dir(simulator)
        commands = self.get_commands()

        task = super().simulate(
            simulator=simulator,
            machine_group=machine_group,
            storage_dir=storage_dir,
            commands=commands,
            sim_config_filename=self.get_config_filename(simulator))

        return task

    @singledispatchmethod
    def set_template_dir(self, simulator: simulators.Simulator):
        pass

    @singledispatchmethod
    def get_config_filename(self, simulator: simulators.Simulator):
        pass

    @singledispatchmethod
    def config_params(self, simulator: simulators.Simulator, input_dir):
        pass

    @singledispatchmethod
    def add_extra_input_files(self, simulator: simulators.Simulator, input_dir):
        pass


@FluidBlock.set_template_dir.register
def _(self, simulator: simulators.SplishSplash):  # pylint: disable=unused-argument
    """Set the template directory for DualSPHysics."""

    self.template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                           SPLISHPLASH_TEMPLATE_INPUT_DIR)


@FluidBlock.get_config_filename.register
def _(cls, simulator: simulators.SplishSplash):  # pylint: disable=unused-argument
    """Returns the configuration filename for SPlisHSPlasH."""
    return SPLISHSPLASH_CONFIG_FILENAME


@FluidBlock.config_params.register
def _(self, simulator: simulators.SplishSplash, input_dir):  # pylint: disable=unused-argument
    """Creates SPlisHSPlasH simulation input files."""

    # Generate the simulation configuration file.
    fluid_margin = 2 * self.params["particle_radius"]

    self.params.update({
        "fluid_position": [
            position + fluid_margin
            for position in self.params["fluid_position"]
        ],
        "fluid_dimensions": [
            dimension - 2 * fluid_margin
            for dimension in self.params["fluid_dimensions"]
        ]
    })


@FluidBlock.add_extra_input_files.register
def _(self, simulator: simulators.SplishSplash, input_dir):  # pylint: disable=unused-argument
    """Add unit box mesh file to input directory."""

    unit_box_file_path = os.path.join(self.template_files_dir, "unit_box.obj")
    shutil.copy(unit_box_file_path, input_dir)


@FluidBlock.set_template_dir.register
def _(self, simulator: simulators.DualSPHysics):  # pylint: disable=unused-argument
    """Set the template directory for DualSPHysics."""

    self.template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                           DUALSPHYSICS_TEMPLATE_INPUT_DIR)


@FluidBlock.get_config_filename.register
def _(cls, simulator: simulators.DualSPHysics):  # pylint: disable=unused-argument
    """Returns the configuration filename for DualSPHysics."""
    return None


@FluidBlock.config_params.register
def _(self, simulator: simulators.DualSPHysics, input_dir):  # pylint: disable=unused-argument
    """Creates DualSPHysics simulation input files."""

    self.params["particle_distance"] = 2 * self.params["particle_radius"]
