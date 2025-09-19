"""Elmer module of the API for multiphysics simulations."""
from typing import List, Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command


@simulators.simulator.mpi_enabled
class Elmer(simulators.Simulator):
    """Class to invoke a generic Elmer simulation on the API.
    """

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the Elmer simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "elmer"

    def run(self,
            input_dir: Optional[str],
            *,
            commands: Optional[List[types.Commands]] = None,
            shell_script: Optional[str] = None,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            on_finish_cleanup: Optional[Union[str, list[str]]] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            commands: List of commands to run using the Elmer simulator.
            shell_script: Path to a shell script (relative to input_dir) to run
                the simulation.
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
            time_to_live: Maximum allowed runtime for the task, specified as a
                string duration. Supports common time duration formats such as
                "10m", "2 hours", "1h30m", or "90s". The task will be
                automatically terminated if it exceeds this duration after
                starting.
            on_finish_cleanup :
                Optional cleanup script or list of shell commands to remove
                temporary or unwanted files generated during the simulation.
                This helps reduce storage usage by discarding unnecessary
                output.
                - If a string is provided, it is treated as the path to a shell
                script that must be included with the simulation files.
                - If a list of strings is provided, each item is treated as an
                individual shell command and will be executed sequentially.
                All cleanup actions are executed in the simulation's working
                directory, after the simulation finishes.
                Examples:
                    on_finish_cleanup = "my_cleanup.sh"

                    on_finish_cleanup = [
                        "rm -rf temp_dir",
                        "rm -f logs/debug.log"
                    ]
            other arguments: See the documentation of the base class.
        """
        if not commands and not shell_script:
            raise ValueError("Either 'commands' or 'shell_script'"
                             " must be provided.")
        if commands and shell_script:
            raise ValueError("Only one of 'commands' or 'shell_script'"
                             " must be provided.")

        if shell_script:

            self._input_files_exist(input_dir=input_dir,
                                    remote_assets=remote_assets,
                                    shell_script=shell_script)
            commands = [f"bash {shell_script}"]

        for i, command in enumerate(commands):
            # Add mpirun if command is a string, contains 'ElmerSolver_mpi', and
            # does not already contain 'mpirun'
            if isinstance(
                    command, str
            ) and "ElmerSolver_mpi" in command and "mpirun" not in command:
                commands[i] = Command(command, mpi_config=on.get_mpi_config())

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           on_finish_cleanup=on_finish_cleanup,
                           **kwargs)
