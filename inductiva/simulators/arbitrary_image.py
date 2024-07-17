"""Class to run commands on an arbitrary image."""
from typing import Optional

from inductiva import types, tasks, simulators


class ArbitraryImage(simulators.Simulator):
    """Class to run commands on an arbitrary image."""

    def __init__(self):
        """Initialize the ArbitraryImage class.
        Point to the API method to run a simulation.
        """
        super().__init__()
        self.api_method_name = "arbitrary.arbitrary_commands.run_simulation"

    def run(self,
            input_dir: str,
            commands: list[str],
            container_image: Optional[str],
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            on: Optional[types.ComputationalResources] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            container_image: The container image to use for the simulation.
                Example: container_image="docker://inductiva/kutu:xbeach_v1.23"
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
                           container_image=container_image,
                           **kwargs)
