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
        input_dir: types.Path,
        sim_config_filename: str,
        output_dir: Optional[types.Path] = None,
        n_cores: int = 1,
        track_logs: bool = False,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           track_logs=track_logs,
                           input_filename=sim_config_filename,
                           n_cores=n_cores,
                           output_dir=output_dir)

    def run_async(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        n_cores: int = 1,
    ) -> str:
        """Run the simulation asynchronously.
        
        Args:
            sim_config_filename: Name of the simulation configuration file.
            n_cores: Number of MPI cores to use for the simulation.
            """

        return super().run_async(input_dir,
                                 n_cores=n_cores,
                                 input_filename=sim_config_filename)
