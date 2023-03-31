"""Describes the physical scenarios and runs its simulation via API."""
from typing import List, Literal, Optional
from enum import Enum
from dataclasses import dataclass

from inductiva.fluids.scenarios.fluid_block import FluidBlock
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.fluid_types import WATER
from inductiva.types import Path

from inductiva.fluids._output_post_processing import SimulationOutput


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

        FluidBlock.__init__(self,
                            density=fluid.density,
                            kinematic_viscosity=fluid.kinematic_viscosity,
                            dimensions=dimensions,
                            position=position)

    #pylint: disable=arguments-renamed
    def simulate(self,
                 device: Literal["cpu", "gpu"] = "cpu",
                 engine: Literal["SPlisHSPlasH",
                                 "DualSPHysics"] = "SPlisHSPlasH",
                 simulation_time: float = 1.,
                 resolution: Literal["high", "medium", "low"] = "medium",
                 output_dir: Optional[Path] = None):
        """Runs SPH simulation of Dam Break scenario.

        Use the FluidBlock class to run a simulation with less freedom of
        choice of the parameters. Namely, the Dam Break scenario only allows
        to change the following arguments.

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
            simulation_time: Simulation time in seconds.
            output_dir: Directory in which the output files will be saved. If
              not specified, the default directory used for API tasks
              (based on an internal ID of the task) will be used.
        """

        particle_radius = ParticleRadius[resolution.upper()].value

        sim_output_path = FluidBlock.simulate(self,
                                              device=device,
                                              engine=engine,
                                              simulation_time=simulation_time,
                                              particle_radius=particle_radius,
                                              output_dir=output_dir)

        return SimulationOutput(sim_output_path)
