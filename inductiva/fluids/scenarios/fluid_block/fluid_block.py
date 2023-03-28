"""Describes the physical scenarios and runs its simulation via API."""
import os
import math
import shutil

from typing import List, Optional

from absl import logging

from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.scenarios.simulator_mixins import SPlisHSPlasHMixin
from inductiva.fluids.scenarios.simulator_mixins import DualSPHysicsMixin
from inductiva.utils.templates import replace_params_in_template

# Global variables to define a scenario
TANK_DIMENSIONS = [1, 1, 1]
TIME_STEP = 0.001
PARTICLE_RADIUS = 0.02
SIMULATION_TIME = 1

SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_block_template.splishsplash.json.jinja"
SPLISHSPLASH_INPUT_FILENAME = "splishsplash_config.json"
UNIT_BOX_MESH_FILENAME = "unit_box.obj"

DUALSPHYSICS_TEMPLATE_FILENAME = "dam_break_template.dualsphysics.xml.jinja"
DUALSPHYSICS_INPUT_FILENAME = "dualsphysics_config.xml"


class FluidBlock(SPlisHSPlasHMixin, DualSPHysicsMixin):
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
            self.position = [0.0, 0.0, 0.0]
        else:
            self.position = position

        if inital_velocity is None:
            self.initial_velocity = [0.0, 0.0, 0.0]
        else:
            self.initial_velocity = inital_velocity

    def _create_aux_files_splishsplash(self, input_dir):
        unit_box_file_path = os.path.join(os.path.dirname(__file__),
                                          UNIT_BOX_MESH_FILENAME)

        shutil.copy(unit_box_file_path, input_dir)

    def _replace_params_in_template_splishsplash(self, input_dir,
                                                 simulator_params):
        """Replaces parameters in the SPlisHSPlasH template file."""

        fluid_margin = 2 * PARTICLE_RADIUS

        replace_params_in_template(
            templates_dir=os.path.dirname(__file__),
            template_filename=SPLISHSPLASH_TEMPLATE_FILENAME,
            params={
                "simulation_time": SIMULATION_TIME,
                "time_step": TIME_STEP,
                "particle_radius": PARTICLE_RADIUS,
                "data_export_rate": 1 / simulator_params.output_time_step,
                "tank_filename": UNIT_BOX_MESH_FILENAME,
                "tank_dimensions": TANK_DIMENSIONS,
                "fluid_filename": UNIT_BOX_MESH_FILENAME,
                "fluid": self.fluid,
                "fluid_position": [fluid_margin] * 3,
                "fluid_dimensions": [
                    dimension - 2 * fluid_margin
                    for dimension in self.dimensions
                ],
                "fluid_velocity": self.initial_velocity,
            },
            output_file_path=os.path.join(input_dir,
                                          SPLISHSPLASH_INPUT_FILENAME),
        )

        logging.info("Estimated number of particles %d",
                     self.estimate_num_particles())
        logging.info("Estimated number of time steps %s",
                     math.ceil(SIMULATION_TIME / TIME_STEP))
        logging.info(
            "Number of output time steps %s",
            math.ceil(SIMULATION_TIME / simulator_params.output_time_step))

    def _replace_params_in_template_dualsphysics(self, input_dir,
                                                 simulator_params):
        """Replaces parameters in the DualSPHysics template file."""

        replace_params_in_template(
            templates_dir=os.path.dirname(__file__),
            template_filename=DUALSPHYSICS_TEMPLATE_FILENAME,
            params={
                "simulation_time": SIMULATION_TIME,
                "particle_distance": 2 * PARTICLE_RADIUS,
                "output_time_step": simulator_params.output_time_step,
                "tank_dimensions": TANK_DIMENSIONS,
                "fluid_dimensions": self.dimensions,
                "fluid_position": self.position,
                "fluid": self.fluid,
            },
            output_file_path=os.path.join(input_dir,
                                          DUALSPHYSICS_INPUT_FILENAME),
        )

    def estimate_num_particles(self):
        """Estimate of the number of SPH particles contained in fluid blocks."""

        # Calculate number of particles for a fluid block
        n_particles_x = round(self.dimensions[0] / (2 * PARTICLE_RADIUS)) - 1
        n_particles_y = round(self.dimensions[1] / (2 * PARTICLE_RADIUS)) - 1
        n_particles_z = round(self.dimensions[2] / (2 * PARTICLE_RADIUS)) - 1

        # Add number of particles to the total sum
        return n_particles_x * n_particles_y * n_particles_z
