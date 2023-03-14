"""Input parameters of SPH simulators."""
from enum import Enum
from dataclasses import dataclass


class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.008
    MEDIUM = 0.012
    LOW = 0.02


@dataclass
class SPHSimulatorParameters:
    output_time_step: float = 1. / 60.
    simulation_time: int = 2


@dataclass
class SPlishSPlasHParameters(SPHSimulatorParameters):
    """Set of parameters for SPLisHSPlasH.

        Args:
            viscosity_solver: Method used to model the viscosity of the fluid.
             The available options are (the default is 'standard'):
                - 'None': Fluid is simulated with no viscosity.
                - 'Standard': Standard SPH formulation of viscosity.
                - 'Weiler-2018': This method is based on the paper "A Physically
            cfl_method: cfl_method: Courant-Friedrichs-Lewy (CFL) method used
              for adaptive time stepping.
              The available options are (default is 'no'):
                - 'no': No adaptive time-stepping is used.
                - 'cfl': Use CFL condition.
                - 'cfl_p': Use CFL condition and consider number of pressure
                  solver iterations."""

    viscosity_solver: str = "Weiler-2018"
    cfl_method: str = "no"


@dataclass
class DualSPHysicsParameters(SPHSimulatorParameters):
    """Set of parameters for DualSPHysics.

        Args:
            cfl_number: Coefficient to multiply dt.
            coefh: Coefficient to calculate the smoothing length in 3D.
            kernel: Interaction Kernel 1:Cubic Spline, 2:Wendland"""
    cflnumber: float = 0.2
    kernel: int = 1
