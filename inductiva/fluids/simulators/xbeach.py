"""DualSPHysics module of the API."""
import pathlib
from typing import Optional

from inductiva.types import Path
from inductiva.fluids.simulators._simulator import Simulator


class XBeach(Simulator):
    """Class to invoke a generic XBeach simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sw.xbeach.run_simulation"

    def simulate(
        self,
        output_dir: Optional[Path] = None,
        n_cores: int = 1,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            device: Device in which to run the simulation.
        """
        return super().simulate(output_dir=output_dir, n_cores=n_cores)
