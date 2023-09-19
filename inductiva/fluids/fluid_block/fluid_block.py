"""Describes the physical scenarios and runs its simulation via API."""

from functools import singledispatchmethod
import os
from typing import List, Optional
import shutil

from inductiva import tasks, resources, fluids, scenarios, simulators, utils

TANK_DIMENSIONS = [1, 1, 1]

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "fluid_block")
SPLISHPLASH_TEMPLATE_INPUT_DIR = "splishsplash"
SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_block_template.splishsplash.json.jinja"
SPLISHSPLASH_CONFIG_FILENAME = "fluid_block.json"
UNIT_BOX_MESH_FILENAME = "unit_box.obj"

DUALSPHYSICS_TEMPLATE_INPUT_DIR = "dualsphysics"
DUALSPHYSICS_TEMPLATE_FILENAME = "fluid_block_template.dualsphysics.xml.jinja"
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

        self.fluid = fluids.FluidType(density=density,
                                      kinematic_viscosity=kinematic_viscosity)

        if len(dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")
        self.dimensions = dimensions

        if position is None:
            self.position = [0, 0, 0]
        else:
            self.position = position

        if initial_velocity is None:
            self.initial_velocity = [0, 0, 0]
        else:
            self.initial_velocity = initial_velocity

    def simulate(
        self,
        simulator: simulators.Simulator = simulators.DualSPHysics(),
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        particle_radius: float = 0.02,
        simulation_time: float = 1,
        adaptive_time_step: bool = True,
        particle_sorting: bool = True,
        time_step: float = 0.001,
        output_time_step: float = 1 / 60,
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
            output_time_step: Time step between outputs, in seconds.
            run_async: Whether to run the simulation asynchronously.
        """
        simulator.override_api_method_prefix("fluid_block")

        # TODO: Avoid storing these as class attributes.
        self.particle_radius = particle_radius
        self.simulation_time = simulation_time
        self.adaptive_time_step = adaptive_time_step
        self.particle_sorting = particle_sorting
        self.time_step = time_step
        self.output_time_step = output_time_step

        task = super().simulate(
            simulator=simulator,
            machine_group=machine_group,
            run_async=run_async,
            sim_config_filename=self.get_config_filename(simulator))

        return task

    @singledispatchmethod
    def get_config_filename(self, simulator: simulators.Simulator):
        pass

    @singledispatchmethod
    def create_input_files(self, simulator: simulators.Simulator):
        pass


@FluidBlock.get_config_filename.register
def _(cls, simulator: simulators.SplishSplash):  # pylint: disable=unused-argument
    """Returns the configuration filename for SPlisHSPlasH."""
    return SPLISHSPLASH_CONFIG_FILENAME


@FluidBlock.create_input_files.register
def _(self, simulator: simulators.SplishSplash, input_dir):  # pylint: disable=unused-argument
    """Creates SPlisHSPlasH simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      SPLISHPLASH_TEMPLATE_INPUT_DIR)
    # Copy the unit box mesh file to the input directory.
    unit_box_file_path = os.path.join(template_files_dir,
                                      UNIT_BOX_MESH_FILENAME)
    shutil.copy(unit_box_file_path, input_dir)

    # Generate the simulation configuration file.
    fluid_margin = 2 * self.particle_radius

    utils.templates.replace_params(
        template_path=os.path.join(template_files_dir,
                                   SPLISHSPLASH_TEMPLATE_FILENAME),
        params={
            "simulation_time": self.simulation_time,
            "time_step": self.time_step,
            "particle_radius": self.particle_radius,
            "data_export_rate": 1 / self.output_time_step,
            "tank_filename": UNIT_BOX_MESH_FILENAME,
            "tank_dimensions": TANK_DIMENSIONS,
            "fluid_filename": UNIT_BOX_MESH_FILENAME,
            "fluid": self.fluid,
            "fluid_position": [
                position + fluid_margin for position in self.position
            ],
            "fluid_dimensions": [
                dimension - 2 * fluid_margin for dimension in self.dimensions
            ],
            "fluid_velocity": self.initial_velocity,
            "particle_sorting": self.particle_sorting,
            "adaptive_time_step": self.adaptive_time_step,
        },
        output_file=os.path.join(input_dir, SPLISHSPLASH_CONFIG_FILENAME),
    )


@FluidBlock.get_config_filename.register
def _(cls, simulator: simulators.DualSPHysics):  # pylint: disable=unused-argument
    """Returns the configuration filename for DualSPHysics."""
    return DUALSPHYSICS_CONFIG_FILENAME


@FluidBlock.create_input_files.register
def _(self, simulator: simulators.DualSPHysics, input_dir):  # pylint: disable=unused-argument
    """Creates DualSPHysics simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      DUALSPHYSICS_TEMPLATE_INPUT_DIR)
    utils.templates.replace_params(
        template_path=os.path.join(template_files_dir,
                                   DUALSPHYSICS_TEMPLATE_FILENAME),
        params={
            "simulation_time": self.simulation_time,
            "time_step": self.time_step,
            "particle_distance": 2 * self.particle_radius,
            "output_time_step": self.output_time_step,
            "tank_dimensions": TANK_DIMENSIONS,
            "fluid_dimensions": self.dimensions,
            "fluid_position": self.position,
            "fluid": self.fluid,
            "adaptive_time_step": self.adaptive_time_step,
        },
        output_file=os.path.join(input_dir, DUALSPHYSICS_CONFIG_FILENAME),
    )
