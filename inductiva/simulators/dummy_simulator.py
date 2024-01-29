"""Dummy simulator for demonstrating and testing purposes.

This simulator acts as an echo of the arguments passed and can be used to
test the commands passed to the simulator.
"""

from typing import Optional, List

from inductiva import types, tasks, simulators


class DummySimulator(simulators.Simulator):
    """Dummy simulator for demo and testing purposes."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "tester.echo.run_simulation"

    def run(self,
            input_dir: types.Path,
            input_filename: str,
            commands: Optional[List[dict]] = None,
            sleep_time: Optional[float] = 1,
            extra_metadata: Optional[dict] = None) -> tasks.Task:
        """Run a dummy simulation.
        
        Args: 
            input_dir: Path to directory with simulation input files.
            input_filename: Name of the test input file.
            commands: List of commands to run on the simulator.
            sleep_time: Time to sleep before running the commands.
            extra_metadata: Extra metadata to be sent to the backend.
        """

        return super().run(input_dir,
                           input_filename=input_filename,
                           commands=commands,
                           sleep_time=sleep_time,
                           extra_metadata=extra_metadata)
