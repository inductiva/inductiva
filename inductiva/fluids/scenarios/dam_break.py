"""Describes the physical scenarios and runs its simulation via API."""
from absl import logging
from typing import List, Literal, Optional
import os
import math
import tempfile
import numpy as np
from enum import Enum
import xml.etree.ElementTree as ET

import inductiva_sph
from inductiva_sph import sph_core
import inductiva
from inductiva.fluids._fluid_types import WATER
from inductiva.fluids._output_post_processing import SimulationOutput
from inductiva.types import Path

# Global variables to define a scenario
COLUMN_VELOCITY = [0.0, 0.0, 0.0]
OUTPUT_TIME_STEP = 1. / 60.
TANK_DIMENSIONS = [1, 1, 1]
FLUID_DIMENSION_LOWER_BOUNDARY = 0.1
FLUID_DIMENSION_UPPER_BOUNDARY = 1
VISCOSITY_SOLVER = "Weiler-2018"
TIME_MAX = 3
INPUT_XML_PATH = os.path.join(os.path.dirname(__file__), 
                              "xml_files/InputCase.xml")


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
                 fluid_position: Optional[List[float]] = None) -> None:
        """Initializes a `DamBreak` object.

        Args:
            fluid_dimensions: A list containing fluid column dimensions,
              in meters.
            fluid: A fluid type to simulate.
            fluid_position: Position of the fluid column in the tank, in meters.
            """

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

    def simulate(self,
                 engine: Literal["SPlisHSPlasH",
                                 "DualSPHysics"] = "DualSPHysics",
                 resolution: Literal["high", "medium", "low"] = "medium",
                 simulation_time: float = 1.,
                 device: Literal["cpu", "gpu"] = "cpu",
                 cfl_method: Literal["no", "cfl", "cfl_p"] = "no",
                 output_dir: Optional[Path] = None):
        """Runs SPH simulation of the Dam Break scenario.

        Args:
            device: Sets the device for a simulation to be run.
            engine: The software platform to be used for the simulation.
              Available options are (the default is DualSPHysics):
              - SPlisHSPlasH
              - DualSPHysics
            resolution: Sets the fluid resolution to simulate.
              Available options are (the default is "medium"):
              - "high"
              - "medium"
              - "low"
            time_max: Maximum time of simulation, in seconds.
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

        self.particle_radius = ParticleRadius[resolution.upper()].value

        if simulation_time > TIME_MAX:
            raise ValueError(
                f"`simulation_time` cannot exceed {TIME_MAX} seconds.")
        self.simulation_time = simulation_time
        self.cfl_method = cfl_method
        self.device = device
        self.output_dir = output_dir

        # Create a temporary directory to store simulation input files
        input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with

        if engine == "SPlisHSPlasH":
            sim_output_path = self._splishsplash_simulation(
                input_dir=input_temp_dir)
        elif engine == "DualSPHysics":
            sim_output_path = self._dualsphysics_simulation(
                input_dir=input_temp_dir)

        # Delete temporary input directory
        input_temp_dir.cleanup()

        inductiva_sph.splishsplash.io_utils.convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(sim_output_path, "particles"), 
            output_time_step=OUTPUT_TIME_STEP, 
            netcdf_data_dir=os.path.join(sim_output_path, "particles", "netcdf"))
        
        return SimulationOutput(sim_output_path)

    def _splishsplash_simulation(self, input_dir):
        """Runs simulation on SPlisHSPlasH via API.

        Args:
            input_dir: Directory where the input file will be stored.
        """
        # Create a dam break scenario
        scenario = self._create_scenario()

        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=self.simulation_time,
            particle_radius=self.particle_radius,
            cfl_method=self.cfl_method,
            output_time_step=OUTPUT_TIME_STEP,
            viscosity_method=VISCOSITY_SOLVER,
            output_directory=input_dir.name)

        # Create input file
        simulation.create_input_file()
        logging.info("Estimated number of particles %d",
                     self.estimate_num_particles())
        logging.info("Estimated number of time steps %s",
                     math.ceil(self.simulation_time / simulation.time_step))
        logging.info("Number of output time steps %s",
                     math.ceil(self.simulation_time / OUTPUT_TIME_STEP))

        # Invoke API
        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            sim_dir=input_dir.name,
            device=self.device,
            output_dir=self.output_dir)

        return sim_output_path

    def _dualsphysics_simulation(self, input_dir):
        """Runs simulation on DualSPHysics via API.

        Args:
            input_dir: Directory where the input file will be stored.
        """
        # Upload XML file of dam break scenario
        input_file = ET.parse(INPUT_XML_PATH)
        root = input_file.getroot()

        # Set simulation parameters according to user input
        root.find("./execution/parameters/parameter[@key='TimeMax']").set(
            "value", str(self.simulation_time))

        fluid_size = root.find(".//drawbox[boxfill='solid']/size")
        fluid_size.set("x", str(self.fluid_dimensions[0]))
        fluid_size.set("y", str(self.fluid_dimensions[1]))
        fluid_size.set("z", str(self.fluid_dimensions[2]))

        fluid_point = root.find(".//drawbox[boxfill='solid']/point")
        fluid_point.set("x", str(self.fluid_position[0]))
        fluid_point.set("y", str(self.fluid_position[1]))
        fluid_point.set("z", str(self.fluid_position[2]))

        particle_size = root.find(".//definition")
        particle_size.set("dp", str(self.particle_radius * 2))

        # Create input file
        input_file.write(os.path.join(input_dir.name, "InputCase.xml"))

        # Invoke API
        return inductiva.sph.dualsphysics.run_simulation(
            sim_dir=input_dir.name,
            input_filename="InputCase",
            device=self.device,
            output_dir=self.output_dir)

    def _create_scenario(self):
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
