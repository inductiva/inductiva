"""Dummy simulator for demonstrating and testing purposes.

This simulator acts as an echo of the arguments passed and can be used to
test the commands passed to the simulator.
"""

from typing import Optional

from inductiva import types, tasks, simulators


class DummySimulator(simulators.Simulator):
    """Dummy simulator for demo and testing purposes."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "tester.echo.run_simulation"

    def run(self,
            input_dir: types.Path,
            input_filename: str,
            sleep_time: Optional[float] = 1,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            extra_metadata: Optional[dict] = None,
            **kwargs) -> tasks.Task:
        """Run a dummy simulation that echo's to a file.
        
        Args: 
            input_dir: Path to directory with simulation input files.
            input_filename: Name of the test input file.
            sleep_time: Time to sleep before running the commands.
            extra_metadata: Extra metadata to be sent to the backend.
        """

        return super().run(input_dir,
                           on=on,
                           input_filename=input_filename,
                           sleep_time=sleep_time,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           **kwargs)
