"""Class to run SeisSol simulations on a custom container image."""

from typing import List, Optional
from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class SeisSol(simulators.Simulator):
    """Class to run SeisSol simulations on a custom container image."""

    def __init__(self, container_image: str):
        """
        Initialize the SeisSol simulator class.

        Args:
            container_image: The container image to use for the simulation.
                Example: container_image="docker://inductiva/seissol:latest"
        """
        self.container_image = container_image
        super().__init__()
        self.simulator = "seissol"
        
    def _get_image_uri(self):
        """Return the URI of the container image."""
        return self.container_image

    def run(self,
            input_dir: Optional[str],
            commands: List[str],
            *,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            **kwargs) -> tasks.Task:
        """
        Run the SeisSol simulation.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run (e.g., SeisSol binaries with parameters).
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
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           container_image=self._image_uri,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
