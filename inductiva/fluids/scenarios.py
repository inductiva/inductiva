"""Describes the physical scenarios and runs its simulation via API."""
import inductiva_sph
from inductiva_sph import sph_core

# Glabal variables to define a scenario
TIME_MAX = 0.6
PARTICLE_RADIUS = 0.002
COLUMN_VELOCITY = [0.0, 0.0, 0.0]
BOUNDARY_RESOLUTION = [20, 20, 20]
OUTPUT_TIME_STEP = 1. / 60.
TANK_LENGTH = 0.84
TANK_WIDTH = 0.5715
TANK_HEIGHT = 0.12
TANK_DIMENSION = [0.84, 0.05715, 0.12]
COLUMN_POSITION = [0.0, 0.0, 0.0]


class DamBreak:
    """Physical scenario of a dam break simulation."""

    def __init__(self, fluid: sph_core.fluids.FluidProperties,
                 fluid_dimensions: list[float, float, float]) -> None:
        """Initializes a `DamBreak` object.

        Args:
            fluid: A fluid type of the simulation.
            fluid_dimensions: A list containing the fluid column dimensions
            relative to the tank dimensions."""

        self.fluid = fluid

        #  Set fluid block dimensions according to the input
        if max(fluid_dimensions) <= 1:
            self.fluid_dimension = [
                fluid_dimensions[0] * TANK_LENGTH,
                fluid_dimensions[1] * TANK_WIDTH,
                fluid_dimensions[2] * TANK_HEIGHT
            ]

    def simulate(self):
        """Runs SPH simulation of the Dam Break simulation."""

        # Create a dam break scenario
        scenario = self.__create_scenario()

        # Create simulation
        # pylint: disable=unused-variable
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=TIME_MAX,
            particle_radius=PARTICLE_RADIUS,
            boundary_resolution=BOUNDARY_RESOLUTION,
            output_time_step=OUTPUT_TIME_STEP)

        #  TODO
        # Create input JSON
        # input_json = simulation.create_input_json()
        # Invoke API
        # inductiva.sph(input_json)

    def __create_scenario(self):

        # Create fluid column
        fluid_block = sph_core.fluids.BoxFluidBlock(
            fluid_properties=self.fluid,
            position=COLUMN_POSITION,
            dimensions=self.fluid_dimension,
            initial_velocity=COLUMN_VELOCITY)

        # Set up scenario
        scenario = sph_core.scenarios.DamBreakSPHScenario(
            dimensions=TANK_DIMENSION, fluid_blocks=[fluid_block])

        return scenario
