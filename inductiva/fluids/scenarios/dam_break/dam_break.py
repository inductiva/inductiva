"""Describes the physical scenarios and runs its simulation via API."""
from typing import List, Literal, Optional
from enum import Enum
from dataclasses import dataclass
from uuid import UUID

from inductiva.simulation import Simulator
from inductiva.fluids.simulators import DualSPHysics
from inductiva.fluids.scenarios.fluid_block import FluidBlock
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.fluid_types import WATER
from inductiva.types import Path

from inductiva.fluids.scenarios._post_processing import SPHSimulationOutput


@dataclass
class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.008
    MEDIUM = 0.012
    LOW = 0.02


class DamBreak(FluidBlock):
    """Physical scenario of a dam break simulation."""

    def __init__(
        self,
        fluid: FluidType = WATER,
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
            dimensions = [0.3, 1, 1]

        super().__init__(density=fluid.density,
                         kinematic_viscosity=fluid.kinematic_viscosity,
                         dimensions=dimensions,
                         position=position)

    # pylint: disable=arguments-renamed
    def simulate(
        self,
        simulator: Simulator = DualSPHysics(),
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        device: Literal["cpu", "gpu"] = "gpu",
        resolution: Literal["high", "medium", "low"] = "medium",
        simulation_time: float = 1,
    ):
        """Simulates the scenario.

        Args:
            simulator: Simulator to use.
            output_dir: Directory to store the simulation output.
            device: Device in which to run the simulation.
            resolution: Resolution of the simulation.
            simulation_time: Simulation time, in seconds.
        """

        particle_radius = ParticleRadius[resolution.upper()].value

        sim_output_path = super().simulate(simulator=simulator,
                                           output_dir=output_dir,
                                           resource_pool_id=resource_pool_id,
                                           device=device,
                                           particle_radius=particle_radius,
                                           simulation_time=simulation_time)

        return SPHSimulationOutput(sim_output_path.sim_output_dir)
