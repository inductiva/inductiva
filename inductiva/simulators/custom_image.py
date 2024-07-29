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

        self.api_method_name = "arbitrary.arbitrary_commands.run_simulation"

    def _get_image_uri(self):
        return self.container_image

    def run(self,
            input_dir: str,
            commands: List[str],
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            on: Optional[types.ComputationalResources] = None,
            resubmit_on_preemption: bool = False,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Parent directory for storing simulation
                               results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           container_image=self._image_uri,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
