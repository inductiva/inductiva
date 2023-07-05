"""Describes the physical scenarios and runs its simulation via API."""

from functools import singledispatchmethod
import os
from typing import List, Literal, Optional
import shutil
from uuid import UUID

from inductiva.types import Path
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.simulators import SPlisHSPlasH
from inductiva.fluids.simulators import DualSPHysics
from inductiva.utils.templates import replace_params_in_template
from inductiva.fluids.scenarios._post_processing import SPHSimulationOutput

TANK_DIMENSIONS = [1, 1, 1]

SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_block_template.splishsplash.json.jinja"
SPLISHSPLASH_CONFIG_FILENAME = "fluid_block.json"
UNIT_BOX_MESH_FILENAME = "unit_box.obj"

DUALSPHYSICS_TEMPLATE_FILENAME = "fluid_block_template.dualsphysics.xml.jinja"
DUALSPHYSICS_CONFIG_FILENAME = "dam_break.xml"


class FluidBlock(Scenario):
    """Physical scenario of a general fluid block simulation."""

    def __init__(self,
                 density: float,
                 kinematic_viscosity: float,
                 dimensions: List[float],
                 position: Optional[List[float]] = None,
                 inital_velocity: Optional[List[float]] = None):
        """Initializes a `FluidBlock` object.

        Args:
            density: Density of the fluid in kg/m^3.
            kinematic_viscosity: Kinematic viscosity of the fluid,
                in m^2/s.
            dimensions: A list containing fluid column dimensions,
                in meters.
            position: Position of the fluid column in the tank,
                in meters.
            initial_velocity: Initial velocity of the fluid block
                in the [x, y, z] axes, in m/s.
        """

        self.fluid = FluidType(density=density,
                               kinematic_viscosity=kinematic_viscosity)

        if len(dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")
        self.dimensions = dimensions

        if position is None:
            self.position = [0, 0, 0]
        else:
            self.position = position

        if inital_velocity is None:
            self.initial_velocity = [0, 0, 0]
        else:
            self.initial_velocity = inital_velocity

    def simulate(
        self,
        simulator: Simulator = DualSPHysics(),
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        device: Literal["cpu", "gpu"] = "gpu",
        particle_radius: float = 0.02,
        simulation_time: float = 1,
        adaptive_time_step: bool = True,
        particle_sorting: bool = True,
        time_step: float = 0.001,
        output_time_step: float = 1 / 60,
    ):
        """Simulates the scenario.

        Args:
            simulator: Simulator to use.
            output_dir: Directory to store the simulation output.
            device: Device in which to run the simulation.
            particle_radius: Radius of the fluid particles, in meters.
            simulation_time: Simulation time, in seconds.
            adaptive_time_step: Whether to use adaptive time step.
            particle_sorting: Whether to use particle sorting.
            time_step: Time step, in seconds.
            output_time_step: Time step between outputs, in seconds.
        """

        # TODO: Avoid storing these as class attributes.
        self.particle_radius = particle_radius
        self.simulation_time = simulation_time
        self.adaptive_time_step = adaptive_time_step
        self.particle_sorting = particle_sorting
        self.time_step = time_step
        self.output_time_step = output_time_step

        output_path = super().simulate(
            simulator,
            output_dir=output_dir,
            resource_pool_id=resource_pool_id,
            device=device,
        )

        # TODO: Add any kind of post-processing here, e.g. convert files?
        # convert_vtk_data_dir_to_netcdf(
        #     data_dir=os.path.join(output_path, "vtk"),
        #     output_time_step=SPLISHSPLASH_OUTPUT_TIM_STEP,
        #     netcdf_data_dir=os.path.join(output_path, "netcdf"))

        return SPHSimulationOutput(output_path)

    @singledispatchmethod
    @classmethod
    def get_config_filename(cls, simulator: Simulator):  # pylint: disable=unused-argument
        raise ValueError(
            f"Simulator not supported for `{cls.__name__}` scenario.")

    @singledispatchmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def gen_config(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )


@FluidBlock.get_config_filename.register
def _(cls, simulator: SPlisHSPlasH):  # pylint: disable=unused-argument
    """Returns the configuration filename for SPlisHSPlasH."""
    return SPLISHSPLASH_CONFIG_FILENAME


@FluidBlock.gen_aux_files.register
def _(self, simulator: SPlisHSPlasH, input_dir):  # pylint: disable=unused-argument
    """Generates auxiliary files for SPlisHSPlasH."""
    unit_box_file_path = os.path.join(os.path.dirname(__file__),
                                      UNIT_BOX_MESH_FILENAME)
    shutil.copy(unit_box_file_path, input_dir)


@FluidBlock.gen_config.register
def _(self, simulator: SPlisHSPlasH, input_dir: str):  # pylint: disable=unused-argument
    """Generates the configuration file for SPlisHSPlasH."""

    fluid_margin = 2 * self.particle_radius

    replace_params_in_template(
        templates_dir=os.path.dirname(__file__),
        template_filename=SPLISHSPLASH_TEMPLATE_FILENAME,
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
        output_file_path=os.path.join(input_dir, SPLISHSPLASH_CONFIG_FILENAME),
    )


@FluidBlock.get_config_filename.register
def _(cls, simulator: DualSPHysics):  # pylint: disable=unused-argument
    """Returns the configuration filename for DualSPHysics."""
    return DUALSPHYSICS_CONFIG_FILENAME


@FluidBlock.gen_aux_files.register
def _(self, simulator: DualSPHysics, input_dir):  # pylint: disable=unused-argument
    pass


@FluidBlock.gen_config.register
def _(self, simulator: DualSPHysics, input_dir: str):  # pylint: disable=unused-argument
    """Generates the configuration file for DualSPHysics."""

    replace_params_in_template(
        templates_dir=os.path.dirname(__file__),
        template_filename=DUALSPHYSICS_TEMPLATE_FILENAME,
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
        output_file_path=os.path.join(input_dir, DUALSPHYSICS_CONFIG_FILENAME),
    )
