"""GROMACS module of the API"""

from typing import Optional, List

from inductiva import types, tasks, resources, simulators


class GROMACS(simulators.Simulator):
    """Class to invoke any GROMACS command on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "md.gromacs.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: List[dict],
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run a list of GROMACS commands.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run using the GROMACS simulator.
            machine_group: The machine group to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           commands=commands,
                           run_async=run_async)
