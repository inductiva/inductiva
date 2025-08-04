"""Delft3D module of the API."""
from typing import Optional, Union

from inductiva import simulators, types


class Delft3D(simulators.Simulator):
    """Class to invoke a generic Delft3D simulation on the API."""

    _default_delft3d_device = "cpu"

    def __init__(
        self,
        /,
        version: Optional[str] = None,
        use_dev: bool = False,
    ):
        super().__init__(version, use_dev, self._default_delft3d_device)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "delft3d"

    def run(self,
            input_dir: Optional[str],
            *,
            commands: Optional[types.Commands] = None,
            shell_script: Optional[str] = None,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            on_finish_cleanup: Optional[Union[str, list[str]]] = None,
            **kwargs):
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            storage_dir: Directory for storing results.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
            time_to_live: Maximum allowed runtime for the task, specified as a
                string duration. Supports common time duration formats such as
                "10m", "2 hours", "1h30m", or "90s". The task will be
                automatically terminated if it exceeds this duration after
                starting.
            commands: List of commands to run the simulation. Cannot be used
                with `shell_script`.
            shell_script: Name of the shell script to run the simulation.
                Cannot be used with `commands`.
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
        if (commands is None) == (shell_script is None):
            raise ValueError("You must provide exactly one of 'commands' or "
                             "'shell_script'.")

        if shell_script is not None:
            #Check if the shell script exists

            self._input_files_exist(input_dir=input_dir,
                                    remote_assets=remote_assets,
                                    shell_script=shell_script)
            commands = [f"bash {shell_script}"]

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
