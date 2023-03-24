"""SWASH module of the API."""
from typing import Optional
from inductiva.fluids.simulators._base_simulator import BaseSimulator
from inductiva.types import Path


class SWASH(BaseSimulator):
    """Class to invoke a generic SWASH simulation on the API."""

    def __init__(self, sim_dir: Path, sim_config_filename: str):
        super().__init__(sim_dir, sim_config_filename,
                         "sw.swash.run_simulation")

    def simulate(self, output_dir: Optional[Path] = None, n_cores=1) -> Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
        """
        return super().simulate(n_cores=n_cores, output_dir=output_dir)
