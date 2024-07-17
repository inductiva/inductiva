"""Class to run commands on an custom image."""
from typing import Optional

from inductiva import types, tasks, simulators


class CustomImage(simulators.Simulator):
    """Class to run commands on an custom image."""

    # pylint: disable=W0231
    def __init__(self, container_image: str):
        """Initialize the ArbitraryImage class.
        Point to the API method to run a simulation.
        Args:
            container_image: The container image to use for the simulation.
                Example: container_image="docker://inductiva/kutu:xbeach_v1.23"
        """
        self._image_uri = container_image
        self._version = "N/A"

        self.api_method_name = "arbitrary.arbitrary_commands.run_simulation"

    @property
    def name(self):
        """Get the name of the simulator."""
        return "N/A"

    def run(self,
            input_dir: str,
            commands: list[str],
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            on: Optional[types.ComputationalResources] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Parent directory for storing simulation
                               results.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           container_image=self._image_uri,
                           **kwargs)
