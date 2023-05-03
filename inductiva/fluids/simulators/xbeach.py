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
        input_dir: types.Path,
        sim_config_filename: Optional[str] = "params.txt",
        output_dir: Optional[types.Path] = None,
        n_cores: int = 1,
        track_logs: bool = False,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            n_cores: Number of MPI cores to use for the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            input_filename=sim_config_filename,
            n_cores=n_cores,
            track_logs=track_logs,
            output_dir=output_dir,
        )

    def run_async(
        self,
        input_dir: types.Path,
        sim_config_filename: Optional[str] = "params.txt",
        n_cores: int = 1,
    ) -> str:
        """Run the simulation asynchronously.
        
        Args:
            sim_config_filename: Name of the simulation configuration file.
            n_cores: Number of MPI cores to use for the simulation.
            """

        return super().run_async(input_dir,
                                 input_filename=sim_config_filename,
                                 n_cores=n_cores)
