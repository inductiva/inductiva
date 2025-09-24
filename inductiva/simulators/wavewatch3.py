"""WAVEWATCH III module of the API for numerical simulations."""
from typing import Optional, Union

from inductiva import simulators, tasks, types


@simulators.simulator.mpi_enabled
class WaveWatch3(simulators.Simulator):
    """Class to invoke a generic WAVEWATCH III simulation on the API.
    """

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the WAVEWATCH simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "WAVEWATCH-III"

    def run(self,
            input_dir: Optional[str],
            *,
            commands: types.Commands = None,
            shell_script: Optional[str] = None,
            switch: Optional[str] = None,
            custom_switch: Optional[str] = None,
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
            commands: List of commands to run the simulation. Can be used
                instead of `shell_scriot`.
            shell_script: Path to a shell script that will be executed to run
                the simulation. Can be used instead of `commands`.
            switch: The WAVEWATCH III switch to use for compilation. Some
                examples are: `Ifremer1`, `Ifremer2_pdlib`, `UKMO_reg` etc.
            custom_switch: Custom switch to use for compilation. Must be sent
                with the `input_dir`.
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
        """
        compile_commands = []

        # Check if both commands and shell_script are provided or if neither is
        if (commands is None) == (shell_script is None):
            raise ValueError("You must provide exactly one of 'commands' or "
                             "'shell_script'.")

        # Check if switch and custom_switch are both provided or if neither is
        no_switches = switch is None and custom_switch is None
        both_switches = switch is not None and custom_switch is not None
        if (no_switches) or (both_switches):
            raise ValueError(
                "You must provide either 'switch' or 'custom_switch'.\n"
                "`switch` (Ex: `Ifremer1`, `Ifremer2_pdlib`, `UKMO_reg`) "
                "is a WAVEWATCH III switch and `custom_switch` is "
                "a custom switch provided with the input directory.")

        # Create the compilation commands based on the provided switch
        if switch is not None:
            compile_commands.append(f"bash /home/compile_ww3.sh {switch} false")
        if custom_switch is not None:
            compile_commands.append(
                f"bash /home/compile_ww3.sh {custom_switch} true")

        # Checks if the input files exist
        input_files_kwargs = {
            "input_dir": input_dir,
            "remote_assets": remote_assets,
        }
        if shell_script:
            input_files_kwargs["shell_script"] = shell_script
        if custom_switch:
            input_files_kwargs["custom_switch"] = custom_switch

        self._input_files_exist(**input_files_kwargs)

        commands = commands if commands is not None else [
            f"bash {shell_script}"
        ]
        commands = compile_commands + commands

        return super().run(input_dir,
                           on=on,
                           storage_dir=storage_dir,
                           commands=commands,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           on_finish_cleanup=on_finish_cleanup,
                           **kwargs)
