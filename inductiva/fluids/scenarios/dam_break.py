"""Describes the physical scenarios and runs its simulation via API."""
from enum import Enum
from typing import List, Literal, Optional

import numpy as np

from inductiva_sph import sph_core
from inductiva.fluids.scenarios.fluid_block import FluidBlock
from inductiva.fluids._fluid_types import WATER
from inductiva.types import Path


# Global variables to define a scenario
TANK_DIMENSIONS = [1., 1., 1.]
FLUID_DIMENSION_LOWER_BOUNDARY = 0.1
FLUID_DIMENSION_UPPER_BOUNDARY = 1


class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.008
    MEDIUM = 0.012
    LOW = 0.02


class DamBreak(FluidBlock):
    """Physical scenario of a dam break simulation."""

    def __init__(self,
                 fluid: sph_core.fluids.FluidProperties = WATER,
                 dimensions: Optional[List[float]] = None,
                 position: Optional[List[float]] = None) -> None:
        """Initializes a `DamBreak` object.

        Args:
            fluid: A fluid type to simulate.
            fluid_dimensions: A list containing fluid column dimensions,
              in meters.
            fluid_position: Position of the fluid column in the tank, in meters.
            particle_radius: Radius of the discretization particles, in meters.
              Used to control particle spacing. Smaller particle radius means a
              finer discretization, hence more particles.
            """

        # Initialize a class even if no values are provided
        if position is None:
            position = [0.0, 0.0, 0.0]
        if dimensions is None:
            dimensions = [0.5, 0.5, 0.5]

        if len(dimensions) != 3:
            raise ValueError("`dimensions` must have 3 values.")
        if len(position) != 3:
            raise ValueError("`position` must have 3 values.")

        #  Set fluid block dimensions according to the input
        if max(dimensions) > FLUID_DIMENSION_UPPER_BOUNDARY:
            raise ValueError(f"The values of `fluid_dimensions` cannot exceed \
                {FLUID_DIMENSION_UPPER_BOUNDARY}.")
        if min(dimensions) < FLUID_DIMENSION_LOWER_BOUNDARY:
            raise ValueError(
                f"The values of `fluid_dimensions` must be larger than \
                {FLUID_DIMENSION_LOWER_BOUNDARY}.")

        if np.greater(np.add(dimensions, position),
                      np.array(TANK_DIMENSIONS)).any():
            raise ValueError("Fluid cannot exceed tank borders.")

        FluidBlock.__init__(self, density=fluid.density,
                            kinematic_viscosity=fluid.kinematic_viscosity,
                            dimensions=dimensions,
                            position=position)

    def simulate(self,
                 device: Literal["cpu", "gpu"] = "cpu",
                 simulation_time: float = 1.,
                 resolution: Literal["high", "medium", "low"] = "medium",
                 output_dir: Optional[Path] = None):
        """Runs SPH simulation of Dam Break scenario.
        
        Use the FluidBlock class to run a simulation with less freedom of
        choice of the parameters. Namely, the Dam Break scenario only allows
        to change the following arguments.

        Args:
            device: Sets the device for a simulation to be run.
            resolution: Sets the fluid resolution to simulate.
            Available options are (the default is "medium"):
            - "high"
            - "medium"
            - "low"
            time_max: Maximum time of simulation, in seconds.
            output_dir: Directory in which the output files will be saved. If
                not specified, the default directory used for API tasks
                (based on an internal ID of the task) will be used.
        """

        particle_radius = ParticleRadius[resolution.upper()].value

        return FluidBlock.simulate(self,
                                   device=device,
                                   particle_radius=particle_radius,
                                   simulation_time=simulation_time,
                                   output_dir=output_dir)
