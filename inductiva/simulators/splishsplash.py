"""SplisHSPlasH simulator module of the API."""
from typing import List, Optional

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
        self.simulator = "splishsplash"

    def run(
        self,
        input_dir: Optional[str],
        sim_config_filename: str,
        *,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[List[str]] = None,
        **kwargs,
    ) -> tasks.Task:
        """Run the SPlisHSPlasH simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Directory for storing simulation results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
        Returns:
            Task object representing the simulation task.
        """
        return super().run(
            input_dir,
            input_filename=sim_config_filename,
            storage_dir=storage_dir,
            on=on,
            resubmit_on_preemption=resubmit_on_preemption,
            remote_assets=remote_assets,
            **kwargs,
        )
