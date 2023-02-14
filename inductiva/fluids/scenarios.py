from inductiva import fluids
import inductiva_sph
from inductiva_sph import sph_core

# Glabal variables to define a scenario
TIME_MAX = 0.6
PARTICLE_RADIUS = 0.002
COLUMN_VELOCITY = [0.0, 0.0, 0.0]
SIMULATION_METHOD = "divergence-free-SPH"
BOUNDARY_HANDLING_METHOD = "volume-maps"
BOUNDARY_RESOLUTION = [20, 20, 20]
OUTPUT_TIME_STEP = 1. / 60.
TANK_LENGTH = 0.84
TANK_WIDTH = 0.5715
TANK_HEIGHT = 0.12
TANK_DIMENTION = [0.84, 0.05715, 0.12]
COLUMN_POSITION = [0.0, 0.0, 0.0]


class DamBreak:
    """Physical scenario of a dam break simulation."""

    def __init__(self, fluid: fluids.WATER, fluid_block: list[float, 
                                                              float]) -> None:
        """Initializes a `DamBreak` object.

        Args:
            fluid: A fluid type of the simulation.
            fluid_blocks: A list containing the fluid block's length and
              height relative to the tank's dimentions."""

        self.fluid = fluid
        self.fluid_block = fluid_block

        #  Set fluid block dimentions according to the input
        if max(fluid_block) <= 1:
            self.column_dimention = [
                fluid_block * TANK_LENGTH, TANK_WIDTH, fluid_block * TANK_HEIGHT
            ]

    def simulate(self):
        """Runs SPH simulation of the Dam Break simulation."""

        # Create a dam break scenario
        scenario = self.__create_scenario()

        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=TIME_MAX,
            particle_radius=PARTICLE_RADIUS,
            simulation_method=SIMULATION_METHOD,
            boundary_handling_method=BOUNDARY_HANDLING_METHOD,
            boundary_resolution=BOUNDARY_RESOLUTION,
            output_time_step=OUTPUT_TIME_STEP)

        # Create input files
        simulation._create_input_file()

        # Invoke API

    def __create_scenario(self):
        fluid_properties = sph_core.fluids.FluidProperties(
            density=fluids.WATER.density,
            kinematic_viscosity=fluids.WATER.kinematic_viscosity)

        # Create fluid column
        fluid_block = sph_core.fluids.BoxFluidBlock(
            fluid_properties=fluid_properties,
            position=COLUMN_POSITION,
            dimensions=self.column_dimention,
            initial_velocity=COLUMN_VELOCITY)

        # Set up scenario
        scenario = sph_core.scenarios.DamBreakSPHScenario(
            dimensions=TANK_DIMENTION, fluid_blocks=[fluid_block])

        return scenario
