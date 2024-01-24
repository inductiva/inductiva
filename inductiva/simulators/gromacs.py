"""GROMACS module of the API"""

from typing import Optional

from inductiva import types, tasks, simulators


class GROMACS(simulators.Simulator):
    """Class to invoke any GROMACS command on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "md.gromacs.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: types.Commands,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        extra_metadata: Optional[dict] = None,
    ) -> tasks.Task:
        """Run a list of GROMACS commands.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run using the GROMACS simulator.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Parent directory for storing simulation
                               results.
        """

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata)
