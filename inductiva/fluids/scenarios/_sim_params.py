"""Input parameters of SPH simulators."""
from enum import Enum
from dataclasses import dataclass


class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.008
    MEDIUM = 0.012
    LOW = 0.02


@dataclass
class SPlishSPlasHParameters:
    """Set of parameters for SPLisHSPlasH.

        Args:
            viscosity_solver: Method used to model the viscosity of the fluid.
             The available options are (the default is 'standard'):
                - 'None'
                - 'Standard'
                - 'Weiler-2018'
            output_time_step: Data is exported and saved every
              'output_time_step' seconds.
            cfl_method: cfl_method: Courant-Friedrichs-Lewy (CFL) method used
              for adaptive time stepping.
              The available options are (default is 'no'):
                - 'no'
                - 'cfl'
                - 'cfl_p'
            simulation_method: Pressure solver to use.
            boundry_handling_method: Method to handle boundary interactions
              with particles. The available options are:
                - 'particle-based'
                - 'volume-maps' """

    viscosity_solver: str = "Weiler-2018"
    output_time_step: float = 1. / 60.
    cfl_method: str = "no"
    simulation_method: str = "divergence-free-SPH"
    boundry_handling_method: str = "particle-based"


@dataclass
class DualSPHysicsParameters:
    """Set of parameters for DualSPHysics.

        Args:
            cfl_number: Coefficient to multiply dt.
            time_out: Time step to export the data."""
    cflnumber: float = 0.2
    time_out: float = 1. / 60.
