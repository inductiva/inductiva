"""SplisHSPlasH simulator module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


class SplishSplash(simulators.Simulator):
    """Class to invoke a generic SPlisHSPlasH simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SPlisHSplasH simulator.
        
        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = "sph.splishsplash.run_simulation"

    def run(
        self,
        input_dir: types.PathOrStr,
        sim_config_filename: str,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.PathOrStr] = "",
        extra_metadata: Optional[dict] = None,
        **kwargs,
    ) -> tasks.Task:
        """Run the SPlisHSPlasH simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Directory for storing simulation results.
        Returns:
            Task object representing the simulation task.
        """
        return super().run(
            input_dir,
            input_filename=sim_config_filename,
            extra_metadata=extra_metadata,
            storage_dir=storage_dir,
            on=on,
            **kwargs,
        )
