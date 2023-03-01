"""Describes the physical scenarios and runs its simulation via API."""
from absl import logging
import tempfile
import numpy as np
from enum import Enum
import math

from typing import List, Optional, Literal

import inductiva
import inductiva_sph
from inductiva_sph import sph_core
from inductiva.fluids._output_post_processing import SimulationOutput
from inductiva.fluids._fluid_types import WATER
from inductiva.types import Path

# Glabal variables to define a scenario
COLUMN_VELOCITY = [0.0, 0.0, 0.0]
OUTPUT_TIME_STEP = 1. / 60.
TANK_DIMENSIONS = [1, 1, 1]
FLUID_DIMENSION_LOWER_BOUNDARY = 0.1
FLUID_DIMENSION_UPPER_BOUNDARY = 1
VISCOSITY_SOLVER = "Weiler-2018"
TIME_MAX = 5


class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.008
    MEDIUM = 0.012
    LOW = 0.02


class DamBreak:
    """Physical scenario of a dam break simulation."""

    def __init__(self,
                 fluid_dimensions: List[float],
                 fluid: sph_core.fluids.FluidProperties = WATER,
                 fluid_position: Optional[List[float]] = None,
                 resolution: Literal["high", "medium", "low"] = "medium",
                 simulation_time: float = 1) -> None:
        """Initializes a `DamBreak` object.

        Args:
            fluid_dimensions: A list containing fluid column dimensions,
              in meters.
            fluid: A fluid type to simulate.
            fluid_position: Position of the fluid column in the tank, in meters.
            particle_radius: Radius of the discretization particles, in meters.
              Used to control particle spacing. Smaller particle radius means a
              finer discretization, hence more particles.
            resolution: Sets the fluid resolution to simulate.
              Available options are (the default is "medium"):
              - "high"
              - "medium"
              - "low"
            simulation_time: Simulation time in seconds."""

        self.fluid = fluid

        #  Set fluid block dimensions according to the input
        if max(fluid_dimensions) > FLUID_DIMENSION_UPPER_BOUNDARY:
            raise ValueError(f"The values of `fluid_dimensions` cannot exceed \
                {FLUID_DIMENSION_UPPER_BOUNDARY}.")
        if min(fluid_dimensions) < FLUID_DIMENSION_LOWER_BOUNDARY:
            raise ValueError(
                f"The values of `fluid_dimensions` must be larger than \
                {FLUID_DIMENSION_LOWER_BOUNDARY}.")
        if len(fluid_dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")

        self.fluid_dimensions = fluid_dimensions

        if fluid_position is None:
            self.fluid_position = [0.0, 0.0, 0.0]

        if len(fluid_position) != 3:
            raise ValueError("`fluid_position` must have 3 values.")

        if np.greater(np.add(self.fluid_dimensions, fluid_position),
                      np.array(TANK_DIMENSIONS)).any():
            raise ValueError("Fluid cannot exceed tank borders.")
        self.fluid_position = fluid_position

        self.particle_radius = ParticleRadius[resolution.upper()].value

        if simulation_time > TIME_MAX:
            raise ValueError(f"`simulation_time` cannot exceed {TIME_MAX} seconds.")
        self.simulation_time = simulation_time

    def simulate(self, output_dir: Optional[Path] = None):
        """Runs SPH simulation of the Dam Break scenario.

        Args:
            output_dir: Directory in which the output files will be saved. If
                not specified, then the default directory used for API tasks
                (based on an internal ID of the task) will be used.
        """

        # Create a dam break scenario
        scenario = self.__create_scenario()

        # Create a temporary directory to store simulation input files
        input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with
        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=self.simulation_time,
            particle_radius=self.particle_radius,
            output_time_step=OUTPUT_TIME_STEP,
            viscosity_method=VISCOSITY_SOLVER,
            output_directory=input_temp_dir.name)

        # Create input file
        simulation.create_input_file()
        logging.info("Estimated number of particles %d",
                     self.estimate_num_particles())
        logging.info("Estimated number of time steps %s", 
                     math.ceil(self.simulation_time /simulation.time_step))
        logging.info("Number of output time steps %s",
                     math.ceil(self.simulation_time / OUTPUT_TIME_STEP))

        logging.info("Running SPlisHSPlasH simulation.")
        # Invoke API
        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            input_temp_dir.name, output_dir=output_dir)
        simulation._output_directory = sim_output_path  #pylint: disable=protected-access

        simulation._convert_output_files(False)  #pylint: disable=protected-access

        # Delete temporary input directory
        input_temp_dir.cleanup()

        return SimulationOutput(sim_output_path)

    def __create_scenario(self):

        # Create fluid column
        fluid_block = sph_core.fluids.BoxFluidBlock(
            fluid_properties=self.fluid,
            position=self.fluid_position,
            dimensions=self.fluid_dimensions,
            initial_velocity=COLUMN_VELOCITY)

        # Set up scenario
        scenario = sph_core.scenarios.DamBreakSPHScenario(
            dimensions=TANK_DIMENSIONS, fluid_blocks=[fluid_block])

        return scenario

    def estimate_num_particles(self):
        """Estimate of the number of SPH particles contained in fluid blocks."""

        # Calculate number of particles for a fluid block
        n_particles_x = round(self.fluid_dimensions[0] /
                              (2 * self.particle_radius)) - 1
        n_particles_y = round(self.fluid_dimensions[1] /
                              (2 * self.particle_radius)) - 1
        n_particles_z = round(self.fluid_dimensions[2] /
                              (2 * self.particle_radius)) - 1

        # Add number of particles to the total sum
        return n_particles_x * n_particles_y * n_particles_z
