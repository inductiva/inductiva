"""DualSPHysics module of the API."""
import pathlib
from typing import Literal, Optional
from dataclasses import dataclass

from inductiva.types import Path
from inductiva.fluids.simulators._simulator import Simulator


class SPlisHSPlasH(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sph.splishsplash.run_simulation"

    def simulate(
        self,
        output_dir: Optional[Path] = None,
        device: Literal["gpu", "cpu"] = "cpu",
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            device: Device in which to run the simulation.
        """
        return super().simulate(output_dir=output_dir, device=device)


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
