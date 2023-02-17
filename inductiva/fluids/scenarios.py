"""Describes the physical scenarios and runs its simulation via API."""
import tempfile
import numpy as np

import inductiva
from inductiva.types import DirPath
from ._output_post_processing import SimulationOutput
import inductiva_sph
from inductiva_sph import sph_core

# Glabal variables to define a scenario
TIME_MAX = 0.6
COLUMN_VELOCITY = [0.0, 0.0, 0.0]
OUTPUT_TIME_STEP = 1. / 60.
TANK_LENGTH = 1
TANK_WIDTH = 1
TANK_HEIGHT = 1
TANK_DIMENSION = [TANK_LENGTH, TANK_WIDTH, TANK_HEIGHT]


class DamBreak:
    """Physical scenario of a dam break simulation."""

    def __init__(self, fluid: sph_core.fluids.FluidProperties,
                 fluid_dimensions: list[float, float, float],
                 fluid_position: list[float, float,
                                      float], particle_radius: float) -> None:
        """Initializes a `DamBreak` object.

        Args:
            fluid: A fluid type of the simulation. Ex.: fluids.WATER
            fluid_dimensions: A list containing the fluid column dimensions
            relative to the tank dimensions.
            fluid_position: Position of the fluid column in the tank.
            particle_radius: Radius of the discretization particles, in meters.
              Used to control particle spacing. Smaller particle radius means a
              finer discretization, hence more particles. """

        self.fluid = fluid

        #  Set fluid block dimensions according to the input
        if max(fluid_dimensions) > 1:
            raise ValueError(
                "The values of `fluid_dimensions` cannot exceed 1.")
        if len(fluid_dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")

        self.fluid_dimensions = [
            fluid_dimensions[0] * TANK_LENGTH, fluid_dimensions[1] * TANK_WIDTH,
            fluid_dimensions[2] * TANK_HEIGHT
        ]

        if particle_radius < 0.01:
            raise ValueError("`particle_radius` must be larger than 0.01.")
        self.particle_radius = particle_radius

        if len(fluid_position) != 3:
            raise ValueError("`fluid_position` must have 3 values.")
        if np.greater(np.add(self.fluid_dimensions, fluid_position),
                      np.array(TANK_DIMENSION)).any():
            raise ValueError("Fluid cannot exceed tank borders.")
        self.fluid_position = fluid_position

    def simulate(self):
        """Runs SPH simulation of the Dam Break scenario."""

        # Create a dam break scenario
        scenario = self.__create_scenario()

        # Create a temporary directory to store simulation input files
        input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with
        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=TIME_MAX,
            particle_radius=self.particle_radius,
            output_time_step=OUTPUT_TIME_STEP,
            output_directory=input_temp_dir.name)

        # Create input file
        simulation.create_input_file()

        #  Invoke API
        sim_output_path = inductiva.sph.run_simulation(
            DirPath(input_temp_dir.name))
        simulation._output_directory = sim_output_path.path  #pylint: disable=protected-access

        simulation._convert_output_files()  #pylint: disable=protected-access

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
            dimensions=TANK_DIMENSION, fluid_blocks=[fluid_block])

        return scenario
