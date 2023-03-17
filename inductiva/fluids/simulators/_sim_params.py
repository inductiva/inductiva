"""Input parameters of SPH simulators."""
from dataclasses import dataclass


@dataclass
class SPlisHSPlasHParameters:
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
        z_sort: Enable z-sort, i.e. periodic particle sorting according to
              their z position. Improves cache hits and therefore the
              performance of the simulation.
        simulation_method: Pressure solver to use.
        boundary_handling_method: Method to handle boundary interactions
          with particles. The available options are:
            - 'particle-based'
            - 'volume-maps'
    """

    viscosity_solver: str = "Weiler-2018"
    output_time_step: float = 1. / 60.
    cfl_method: str = "no"
    z_sort: bool = False
    simulation_method: str = "divergence-free-SPH"
    boundary_handling_method: str = "particle-based"


@dataclass
class DualSPHysicsParameters:
    """Set of parameters for DualSPHysics.

    Args:
        cfl_number: Coefficient to multiply dt.
        time_out: Time step to export the data.
    """
    cflnumber: float = 0.2
    output_time_step: float = 0.01
