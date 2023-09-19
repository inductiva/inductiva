"""Describes the physical scenarios and runs its simulation via API."""
from typing import List, Literal, Optional
from enum import Enum
from dataclasses import dataclass

from inductiva import tasks, resources, simulators, fluids


@dataclass
class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.008
    MEDIUM = 0.012
    LOW = 0.02


class DamBreak(fluids.FluidBlock):
    """Physical scenario of a dam break simulation."""

    def __init__(
        self,
        fluid: fluids.FluidType = fluids.WATER,
        dimensions: Optional[List[float]] = None,
        position: Optional[List[float]] = None,
    ):
        """Initializes a `DamBreak` object.

        Args:
            fluid: A fluid type to simulate.
            dimensions: A list containing fluid column dimensions,
              in meters.
            pisition: A list containing fluid column position in a tank,
              in meters.
            """

        if dimensions is None:
            dimensions = [0.3, 0.3, 0.3]

        super().__init__(density=fluid.density,
                         kinematic_viscosity=fluid.kinematic_viscosity,
                         dimensions=dimensions,
                         position=position)

    # pylint: disable=arguments-renamed
    def simulate(
        self,
        simulator: simulators.Simulator = simulators.DualSPHysics(),
        machine_group: Optional[resources.MachineGroup] = None,
        resolution: Literal["high", "medium", "low"] = "low",
        simulation_time: float = 1,
        run_async: bool = False,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: Simulator to use.
            machine_group: The machine group to use for the simulation.
            resolution: Resolution of the simulation.
            simulation_time: Simulation time, in seconds.
            run_async: Whether to run the simulation asynchronously.
        """
        simulator.override_api_method_prefix("dam_break")

        self.particle_radius = ParticleRadius[resolution.upper()].value
        self.simulation_time = simulation_time
        self.adaptive_time_step = True
        self.particle_sorting = True
        self.time_step = 0.001
        self.output_time_step = 1 / 60

        # Inherit the simulate from the Parent of FluidBlock (Scenario) to
        # avoid overriding the api_method_prefix with the one of FluidBlock.
        task = super(fluids.FluidBlock, self).simulate(
            simulator=simulator,
            machine_group=machine_group,
            run_async=run_async,
            sim_config_filename=self.get_config_filename(simulator))

        return task
