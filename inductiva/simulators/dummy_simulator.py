"""Dummy simulator for demonstrating and testing purposes.

This simulator acts as an echo of the arguments passed and can be used to
test the commands passed to the simulator.
"""

from typing import Optional

from inductiva import types, tasks, simulators


class DummySimulator(simulators.Simulator):
    """Dummy simulator for demo and testing purposes."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = "tester.echo.run_simulation"

    def _get_image_uri(self):
        return None

    def run(self,
            input_dir: str,
            input_filename: str,
            sleep_time: Optional[float] = 1,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            extra_metadata: Optional[dict] = None,
            **kwargs) -> tasks.Task:
        """Run a dummy simulation that echo's to a file.

        Args:
            input_dir: Path to directory with simulation input files.
            input_filename: Name of the test input file.
            sleep_time: Time to sleep before running the commands.
            extra_metadata: Extra metadata to be sent to the backend.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
        """

        return super().run(input_dir,
                           on=on,
                           input_filename=input_filename,
                           sleep_time=sleep_time,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
