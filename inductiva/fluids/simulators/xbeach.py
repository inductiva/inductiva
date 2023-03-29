"""DualSPHysics module of the API."""
import pathlib
from typing import Optional

from inductiva import types
from inductiva.simulation import Simulator


class XBeach(Simulator):
    """Class to invoke a generic XBeach simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "sw.xbeach.run_simulation"

    def run(
        self,
        sim_dir: types.Path,
        output_dir: Optional[types.Path] = None,
        n_cores: int = 1,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
        """
        return super().run(sim_dir, output_dir=output_dir, n_cores=n_cores)
