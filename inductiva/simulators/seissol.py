"""Class to run SeisSol simulations on a custom container image."""

from typing import List, Optional
from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class SeisSol(simulators.Simulator):
    """Class to run SeisSol simulations on a custom container image."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the CaNS simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "seissol"

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            n_vcpus: Optional[int] = None,
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            **kwargs) -> tasks.Task:
        """
        Run the SeisSol simulation.

        Args:
            input_dir: Path to the directory containing the input files.
            on: The computational resource to launch the simulation on.
            storage_dir: Parent directory for storing simulation results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
        """
        return super().run(input_dir,
                           n_vcpus=n_vcpus,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
