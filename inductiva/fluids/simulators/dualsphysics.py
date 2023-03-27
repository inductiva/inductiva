"""DualSPHysics module of the API."""
from dataclasses import dataclass
from typing import Literal

import inductiva
from inductiva.types import Path
from inductiva.fluids.simulators._base_simulator import BaseSimulator


class DualSPHysics(BaseSimulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def simulate(self,
                 device: Literal["gpu", "cpu"] = "cpu",
                 output_dir=None) -> Path:
        """Run the simulation.

        Args:
            output_dir: Directory where the generated files will be stored.
        """
        return inductiva.sph.dualsphysics.run_simulation(self.sim_dir,
                                                         self.input_filename,
                                                         device,
                                                         output_dir=output_dir)


@dataclass
class DualSPHysicsParameters:
    """Set of parameters for DualSPHysics.

    Args:
        cfl_number: Coefficient to multiply dt.
        time_out: Time step to export the data.
    """
    cflnumber: float = 0.2
    output_time_step: float = 0.01
