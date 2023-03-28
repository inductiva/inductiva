"""SWASH module of the API."""
from typing import Optional
from inductiva.fluids.simulators._simulator import Simulator
from inductiva.types import Path


class SWASH(Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sw.swash.run_simulation"

    def simulate(self, output_dir: Optional[Path] = None, n_cores=1) -> Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
        """
        return super().simulate(n_cores=n_cores, output_dir=output_dir)
