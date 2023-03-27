"""Describes the physical scenarios and runs its simulation via API."""
import tempfile
import os
import math
from typing import List, Literal, Optional, Union
import shutil
from pathlib import PurePosixPath

from absl import logging

import inductiva
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.simulators import SPlisHSPlasHParameters
from inductiva.fluids.simulators import DualSPHysicsParameters
from inductiva.fluids.post_processing.splishsplash import convert_vtk_data_dir_to_netcdf
from inductiva.types import Path
from inductiva.utils.templates import replace_params_in_template

# Global variables to define a scenario
TANK_DIMENSIONS = [1, 1, 1]
TIME_STEP = 0.001

SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_block_template.splishsplash.json.jinja"
SPLISHSPLASH_INPUT_FILENAME = "fluid_block.json"
UNIT_BOX_MESH_FILENAME = "unit_box.obj"

DUALSPHYSICS_TEMPLATE_FILENAME = "dam_break_template.dualsphysics.xml.jinja"
DUALSPHYSICS_INPUT_FILENAME = "dam_break.xml"


class FluidBlock:
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

    def simulate(self,
                 device: Literal["cpu", "gpu"] = "cpu",
                 engine: Literal["DualSPHysics",
                                 "SPlisHSPlasH"] = "SPlisHSPlasH",
                 simulation_time: float = 1.,
                 particle_radius: float = 0.015,
                 output_dir: Optional[Path] = None,
                 engine_parameters: Union[
                     DualSPHysicsParameters,
                     SPlisHSPlasHParameters] = SPlisHSPlasHParameters):
        """Runs SPH simulation of the fluid block scenario.

        Args:
            device: Sets the device for a simulation to be run.
            engine: The software platform to be used for the simulation.
            Available options are (the default is DualSPHysics):
            - SPlisHSPlasH
            - DualSPHysics
            particle_radius: Radius of the discretization particles, in meters.
            Used to control particle spacing. Smaller particle radius means a
            finer discretization, hence more particles.
            time_max: Maximum time of simulation, in seconds.
            output_dir: Directory in which the output files will be saved. If
                not specified, the default directory used for API tasks
                (based on an internal ID of the task) will be used.
            engine_parameters: Simulator specific parameters.
        """
        self.particle_radius = particle_radius
        self.simulation_time = simulation_time
        self.device = device
        self.output_dir = output_dir
        self.engine_parameters = engine_parameters
        self.fluid_block_dir = os.path.dirname(__file__)

        # Create a temporary directory to store simulation input files
        self.input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with


        if engine.lower() == "splishsplash" and \
            isinstance(engine_parameters, SPlisHSPlasHParameters):
            sim_output_path = self._splishsplash_simulation()
        elif engine.lower() == "dualsphysics" and \
            isinstance(engine_parameters, DualSPHysicsParameters):
            sim_output_path = self._dualsphysics_simulation()
        else:
            raise ValueError("Entered `engine` does not exist or it \
                             does not match with `engine_parameters` class")

        # Delete temporary input directory
        self.input_temp_dir.cleanup()

        return sim_output_path

    def _splishsplash_simulation(self):
        """Runs SPlisHSPlasH simulation via API."""

        input_dir = self.input_temp_dir.name

        unit_box_file_path = os.path.join(self.fluid_block_dir,
                                          UNIT_BOX_MESH_FILENAME)
        shutil.copy(unit_box_file_path, input_dir)

        fluid_margin = 2 * self.particle_radius

        replace_params_in_template(
            templates_dir=self.fluid_block_dir,
            template_filename=SPLISHSPLASH_TEMPLATE_FILENAME,
            params={
                "simulation_time": self.simulation_time,
                "time_step": TIME_STEP,
                "particle_radius": self.particle_radius,
                "data_export_rate": 1 / self.engine_parameters.output_time_step,
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
                     math.ceil(self.simulation_time / TIME_STEP))
        logging.info(
            "Number of output time steps %s",
            math.ceil(self.simulation_time /
                      self.engine_parameters.output_time_step))

        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            sim_dir=input_dir,
            input_filename=SPLISHSPLASH_INPUT_FILENAME,
            device=self.device,
            output_dir=self.output_dir)

        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(sim_output_path, "vtk"),
            output_time_step=self.engine_parameters.output_time_step,
            netcdf_data_dir=os.path.join(sim_output_path, "netcdf"))

        return sim_output_path

    def _dualsphysics_simulation(self):
        """Runs simulation on DualSPHysics via API."""

        replace_params_in_template(
            templates_dir=self.fluid_block_dir,
            template_filename=DUALSPHYSICS_TEMPLATE_FILENAME,
            params={
                "simulation_time": self.simulation_time,
                "particle_radius": 2 * self.particle_radius,
                "output_time_step": self.engine_parameters.output_time_step,
                "tank_dimensions": TANK_DIMENSIONS,
                "fluid_dimensions": self.dimensions,
                "fluid_position": self.position,
                "fluid": self.fluid,
            },
            output_file_path=os.path.join(self.input_temp_dir.name,
                                          DUALSPHYSICS_INPUT_FILENAME),
        )

        return inductiva.sph.dualsphysics.run_simulation(
            sim_dir=self.input_temp_dir.name,
            input_filename=PurePosixPath(DUALSPHYSICS_INPUT_FILENAME).stem,
            device=self.device,
            output_dir=self.output_dir)

    def estimate_num_particles(self):
        """Estimate of the number of SPH particles contained in fluid blocks."""

        # Calculate number of particles for a fluid block
        n_particles_x = round(self.dimensions[0] /
                              (2 * self.particle_radius)) - 1
        n_particles_y = round(self.dimensions[1] /
                              (2 * self.particle_radius)) - 1
        n_particles_z = round(self.dimensions[2] /
                              (2 * self.particle_radius)) - 1

        # Add number of particles to the total sum
        return n_particles_x * n_particles_y * n_particles_z
