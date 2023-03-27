"""DualSPHysics module of the API."""
import pathlib
from typing import Literal, Optional
from dataclasses import dataclass

from inductiva.types import Path
from inductiva.fluids.simulators._simulator import Simulator


class DualSPHysics(Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sph.dualsphysics.run_simulation"

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
class DualSPHysicsParameters:
    """Set of parameters for DualSPHysics.

    Args:
        cfl_number: Coefficient to multiply dt.
        time_out: Time step to export the data.
    """
    cflnumber: float = 0.2
    output_time_step: float = 0.01
