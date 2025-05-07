"""Class to run commands on an custom image."""
from typing import List, Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class CustomImage(simulators.Simulator):
    """Class to run commands on an custom image."""

    def __init__(self, container_image: str):
        """Initialize the ArbitraryImage class.
        Point to the API method to run a simulation.
        Args:
            container_image: The container image to use for the simulation.
                Example: container_image="docker://inductiva/kutu:xbeach_v1.23"
        """
        self.container_image = container_image
        super().__init__()

        self.simulator = "arbitrary_commands"

    def _get_image_uri(self):
        return self.container_image

    def run(self,
            input_dir: Optional[str],
            commands: List[str],
            *,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            project: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            on: The computational resource to launch the simulation on.
            storage_dir: Parent directory for storing simulation
                               results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           container_image=self._image_uri,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
