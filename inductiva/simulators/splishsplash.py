"""SplisHSPlasH simulator module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


class SplishSplash(simulators.Simulator):
    """Class to invoke a generic SPlisHSPlasH simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sph.splishsplash.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        extra_metadata: Optional[dict] = None,
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
            on=on,
            input_filename=sim_config_filename,
            storage_dir=storage_dir,
            extra_metadata=extra_metadata,
        )
