"""Describes the physical scenarios and runs its simulation via API."""
from absl import logging
import tempfile
import math
from typing import List, Literal, Optional
import numpy as np

import inductiva_sph
from inductiva_sph import sph_core
import inductiva
from inductiva.fluids._output_post_processing import SimulationOutput
from inductiva.types import Path

# Global variables to define a scenario
COLUMN_VELOCITY = [0.0, 0.0, 0.0]
OUTPUT_TIME_STEP = 1. / 60.
TANK_DIMENSIONS = [1, 1, 1]
FLUID_DIMENSION_LOWER_BOUNDARY = 0.1
FLUID_DIMENSION_UPPER_BOUNDARY = 1
SIMULATION_METHOD = "divergence-free-SPH"
VISCOSITY_SOLVER = "Weiler-2018"
BOUNDARY_HANDLING_METHOD = "particle-based"


class FluidBlock:
    """Physical scenario of a general fluid block simulation."""

    def __init__(self,
                 fluid_density: float,
                 fluid_kinematic_viscosity: float,
                 fluid_dimensions: List[float],
                 fluid_position: Optional[List[float]]= None,
                 fluid_inital_velocity: Optional[List[float]] = None) -> None:
        """Initializes a `FluidBlock` object.
        
        Args:
            fluid_density: Density of the fluid in kg/m^3.
            fluid_kinematic_viscosity: Kinematic viscosity of the fluid,
                in m^2/s.
            fluid_dimensions: A list containing fluid column dimensions,
                in meters.
            fluid: A fluid type to simulate.
            fluid_position: Position of the fluid column in the tank,
                in meters.
            fluid_initial_velocity: Initial velocity of the fluid block
                in the [x, y, z] axes, in m/s.
        """

        self.fluid = sph_core.fluids.FluidProperties(
            density=fluid_density,
            kinematic_viscosity=fluid_kinematic_viscosity)

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
            fluid_position = [0.0, 0.0, 0.0]

        if np.greater(np.add(self.fluid_dimensions, fluid_position),
                        np.array(TANK_DIMENSIONS)).any():
            raise ValueError("Fluid cannot exceed tank borders.")

        self.fluid_position = fluid_position

        if fluid_inital_velocity is None:
            self.fluid_initial_velocity = [0., 0., 0.]
        else:
            self.fluid_initial_velocity = fluid_inital_velocity

    def simulate(self,
                 device: Literal["cpu", "gpu"] = "cpu",
                 simulation_time: float = 1.,
                 output_time_step: float = 1. / 60.,
                 particle_radius: float = 0.015,
                 cfl_method: Literal["no", "cfl", "cfl_p"] = "no",
                 output_dir: Optional[Path] = None):
        """Runs SPH simulation of the Dam Break scenario.
        Args:
            device: Sets the device for a simulation to be run.
            resolution: Sets the fluid resolution to simulate.
                Available options are (the default is "medium"):
                - "high"
                - "medium"
                - "low"
            time_max: Maximum time of simulation, in seconds.
            Simulation data is exported and saved every
                `output_time_step` seconds.
            cfl_method: cfl_method: Courant-Friedrichs-Lewy (CFL) method used
                for adaptive time stepping. Used to find a time step as large
                as possible to achieve high performance but sufficiently small
                to maintain stability.
                The available options are:
                - 'no': No adaptive time-stepping is used.
                - 'cfl': Use CFL condition.
                - 'cfl_p': Use CFL condition and consider number of pressure
                solver iterations.
            output_dir: Directory in which the output files will be saved. If
                not specified, the default directory used for API tasks
                (based on an internal ID of the task) will be used.
        """

        # Create a dam break scenario
        scenario = self.__create_scenario()

        # Create a temporary directory to store simulation input files
        input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with
        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=simulation_time,
            time_step=output_time_step,
            particle_radius=particle_radius,
            simulation_method=SIMULATION_METHOD,
            viscosity_method = VISCOSITY_SOLVER,
            boundary_handling_method=BOUNDARY_HANDLING_METHOD,
            cfl_method=cfl_method,
            output_time_step=output_time_step,
            output_directory=input_temp_dir.name)

        logging.info(input_temp_dir)

        # Create input file
        simulation.create_input_file()
        logging.info("Estimated number of particles %d",
                     self.estimate_num_particles())
        logging.info("Estimated number of time steps %s",
                     math.ceil(self.simulation_time / simulation.time_step))
        logging.info("Number of output time steps %s",
                     math.ceil(self.simulation_time / OUTPUT_TIME_STEP))

        logging.info("Running SPlisHSPlasH simulation.")
        # Invoke API
        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            sim_dir=input_temp_dir.name, device=device, output_dir=output_dir)
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
            initial_velocity=self.fluid_initial_velocity)

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
