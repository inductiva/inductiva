"""Describes the physical scenarios and runs its simulation via API."""
from typing import List, Literal, Optional, Union
from enum import Enum
from dataclasses import dataclass

from inductiva.fluids.scenarios.fluid_block import FluidBlock
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.fluid_types import WATER
from inductiva.fluids.simulators import SPlisHSPlasHParameters
from inductiva.fluids.simulators import DualSPHysicsParameters
from inductiva.types import Path


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
        super().__init__(density=fluid.density,
                         kinematic_viscosity=fluid.kinematic_viscosity,
                         dimensions=dimensions,
                         position=position)

    # TODO: Remove this method and control *all* parameters via
    # `simulator_params`?
    def simulate(
        self,
        simulator_params: Union[
            DualSPHysicsParameters,
            SPlisHSPlasHParameters] = SPlisHSPlasHParameters(),
        output_dir: Optional[Path] = None,
        resolution: Literal["high", "medium", "low"] = "medium",
    ):
        """Simulates the scenario."""

        # Override the particle radius according to the resolution.
        simulator_params.particle_radius = ParticleRadius[
            resolution.upper()].value

        return super().simulate(simulator_params, output_dir)
