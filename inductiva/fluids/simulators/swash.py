"""SWASH module of the API."""
import pathlib
from typing import Optional
from inductiva.simulation import Simulator
from inductiva import types


class SWASH(Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sw.swash.run_simulation"

    def run(
        self,
        sim_dir: types.Path,
        output_dir: Optional[types.Path] = None,
        n_cores=1,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
        """
        return super().run(sim_dir, n_cores=n_cores, output_dir=output_dir)
