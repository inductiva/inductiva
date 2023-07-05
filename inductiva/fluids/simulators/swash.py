"""SWASH module of the API."""
import pathlib
from typing import Optional
from uuid import UUID
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
        resource_pool_id: Optional[UUID] = None,
        n_cores: int = 1,
    ) -> pathlib.Path:
        """Run the simulation.

        Args:
            n_cores: Number of MPI cores to use for the simulation.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           resource_pool_id=resource_pool_id,
                           input_filename=sim_config_filename,
                           n_cores=n_cores,
                           output_dir=output_dir)

    def run_async(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        resource_pool_id: Optional[UUID] = None,
        n_cores: int = 1,
    ) -> str:
        """Run the simulation asynchronously.

        Args:
            sim_config_filename: Name of the simulation configuration file.
            n_cores: Number of MPI cores to use for the simulation.
            """

        return super().run_async(input_dir,
                                 resource_pool_id=resource_pool_id,
                                 n_cores=n_cores,
                                 input_filename=sim_config_filename)
