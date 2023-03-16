"""Describes the physical scenarios and runs its simulation via API."""
from typing import List, Literal, Optional

from inductiva_sph import sph_core
from inductiva.fluids.scenarios.fluid_block import FluidBlock
from inductiva.fluids._fluid_types import WATER
from inductiva.types import Path
from inductiva.fluids.simulators import ParticleRadius


class DamBreak(FluidBlock):
    """Physical scenario of a dam break simulation."""

    def __init__(self,
                 fluid: sph_core.fluids.FluidProperties = WATER,
                 dimensions: Optional[List[float]] = None):
        """Initializes a `DamBreak` object.

        Args:
            fluid: A fluid type to simulate.
            fluid_dimensions: A list containing fluid column dimensions,
              in meters.
            """
        FluidBlock.__init__(self,
                            density=fluid.density,
                            kinematic_viscosity=fluid.kinematic_viscosity,
                            dimensions=dimensions)

    def simulate(self,
                 engine: Literal["SPlisHSPlasH",
                                 "DualSPHysics"] = "DualSPHysics",
                 resolution: Literal["high", "medium", "low"] = "medium",
                 simulation_time: float = 1.,
                 device: Literal["cpu", "gpu"] = "cpu",
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

        FluidBlock.simulate(self,
                            device=device,
                            engine=engine,
                            simulation_time=simulation_time,
                            particle_radius=particle_radius,
                            output_dir=output_dir)
