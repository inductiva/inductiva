"""Describes the physical scenarios and runs its simulation via API."""
import tempfile
import os
import math
from typing import List, Literal, Optional

from absl import logging
import numpy as np

import inductiva_sph
from inductiva_sph import sph_core
import inductiva
from inductiva.fluids._output_post_processing import SimulationOutput
from inductiva.types import Path

# Global variables to define a scenario
TANK_DIMENSIONS = [1, 1, 1]
SIMULATION_METHOD = "divergence-free-SPH"
VISCOSITY_SOLVER = "Weiler-2018"
BOUNDARY_HANDLING_METHOD = "particle-based"


class FluidBlock:
    """Physical scenario of a general fluid block simulation."""

    def __init__(self,
                 density: float,
                 kinematic_viscosity: float,
                 dimensions: List[float],
                 position: Optional[List[float]]= None,
                 inital_velocity: Optional[List[float]] = None) -> None:
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

        self.fluid = sph_core.fluids.FluidProperties(
            density=density,
            kinematic_viscosity=kinematic_viscosity)

        if len(dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")

        self.dimensions = dimensions

        if position is None:
            self.position = [0.0, 0.0, 0.0]
        else:
            self.position = position

        if np.greater(np.add(self.dimensions, position),
                        np.array(TANK_DIMENSIONS)).any():
            raise ValueError("Fluid cannot exceed tank borders.")

        if inital_velocity is None:
            self.initial_velocity = [0.0, 0.0, 0.0]
        else:
            self.initial_velocity = inital_velocity

    def simulate(self,
                 device: Literal["cpu", "gpu"] = "cpu",
                 simulation_time: float = 1.,
                 output_time_step: float = 1. / 60.,
                 particle_radius: float = 0.015,
                 z_sort: bool = False,
                 cfl_method: Literal["no", "cfl", "cfl_p"] = "no",
                 output_dir: Optional[Path] = None):
        """Runs SPH simulation of the Fluid Block scenario.
        Args:
            device: Sets the device for a simulation to be run.
            particle_radius: Radius of the discretization particles, in meters.
              Used to control particle spacing. Smaller particle radius means a
              finer discretization, hence more particles.
            time_max: Maximum time of simulation, in seconds.
            output_time_step: Simulation data is exported and saved every
                `output_time_step` seconds.
            z_sort: Enable z-sort, i.e. periodic particle sorting according to
              their z position. Improves cache hits and therefore the
              performance of the simulation.
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

        self.particle_radius=particle_radius
        self.simulation_time=simulation_time

        # Create a temporary directory to store simulation input files
        input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with
        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=self.simulation_time,
            time_step=output_time_step,
            particle_radius=self.particle_radius,
            simulation_method=SIMULATION_METHOD,
            viscosity_method=VISCOSITY_SOLVER,
            boundary_handling_method=BOUNDARY_HANDLING_METHOD,
            z_sort=z_sort,
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
                     math.ceil(self.simulation_time / output_time_step))

        logging.info("Running SPlisHSPlasH simulation.")
        # Invoke API
        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            sim_dir=input_temp_dir.name, device=device, output_dir=output_dir)
        simulation._output_directory = sim_output_path  #pylint: disable=protected-access

        # Delete temporary input directory
        input_temp_dir.cleanup()

        inductiva_sph.splishsplash.io_utils.convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(sim_output_path, "vtk"),
            output_time_step=output_time_step,
            netcdf_data_dir=os.path.join(sim_output_path, "netcdf"))

        return SimulationOutput(sim_output_path)

    def __create_scenario(self):
        # Create fluid column
        fluid_block = sph_core.fluids.BoxFluidBlock(
            fluid_properties=self.fluid,
            position=self.position,
            dimensions=self.dimensions,
            initial_velocity=self.initial_velocity)

        # Set up scenario
        scenario = sph_core.scenarios.DamBreakSPHScenario(
            dimensions=TANK_DIMENSIONS, fluid_blocks=[fluid_block])

        return scenario

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
