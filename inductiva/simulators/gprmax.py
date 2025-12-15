"""gprmax module of the API."""
from typing import Optional, Union

from inductiva import simulators, types
from inductiva.commands import MPIConfig
from inductiva.commands.commands import Command


@simulators.simulator.mpi_enabled
class GprMax(simulators.Simulator):
    """Class to invoke a generic gprmax simulation on the API."""

    _default_gprmax_device = "cpu"

    def __init__(
        self,
        /,
        version: Optional[str] = None,
        use_dev: bool = False,
    ):
        super().__init__(version, use_dev, self._default_gprmax_device)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "GprMax"

    def run(self,
            input_dir: Optional[str],
            commands: types.Commands,
            *,
            on: types.ComputationalResources,
            use_gprmax_mpi: bool = False,
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
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
            commands: List of commands to run the simulation.
            use_gprmax_mpi: If True, use GprMax's built-in MPI mechanism
                (requires -mpi flag in command). If False and n_vcpus is
                specified, will use mpirun wrapper (requires --mpi-no-spawn
                flag in command). Default is False.
            use_hwthread: If specified, Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available. Only used when use_gprmax_mpi=False.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                and MPI is being used, all vCPUs will be used.
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

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets)

        # Create MPI config for mpirun wrapper mode
        mpi_config = None
        if not use_gprmax_mpi and n_vcpus is not None:
            mpi_kwargs = {"use_hwthread_cpus": use_hwthread, "np": n_vcpus}
            mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        # Process commands and apply MPI configuration
        # IMPORTANT: MPI is only applied to GprMax simulation commands
        # (those containing "gprMax" but not "tools.").
        # Utility commands like tools.outputfiles_merge or tools.plot_Bscan
        # will NOT have MPI applied, which is the correct behavior since
        # they are post-processing tools that should run single-threaded.
        processed_commands = []
        for idx, command in enumerate(commands):
            if not isinstance(command, str):
                # If it's already a Command object, keep it as is
                processed_commands.append(command)
                continue

            # Check if this is a GprMax simulation command
            # Exclude utility commands (tools.*) from MPI processing
            is_gprmax_cmd = "gprMax" in command and "tools." not in command

            if is_gprmax_cmd:
                # Validate and process GprMax simulation command with MPI
                processed_cmd = self._process_gprmax_command(
                    command, use_gprmax_mpi, mpi_config, n_vcpus, idx)
                processed_commands.append(processed_cmd)
            else:
                # Non-GprMax simulation command (utilities, pre/post-processing)
                # Keep as-is without MPI
                processed_commands.append(command)

        return super().run(input_dir,
                           on=on,
                           commands=processed_commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           on_finish_cleanup=on_finish_cleanup,
                           **kwargs)

    def _process_gprmax_command(self, command: str, use_gprmax_mpi: bool,
                                mpi_config: Optional[MPIConfig],
                                n_vcpus: Optional[int],
                                cmd_idx: int) -> Union[str, Command]:
        """Process a GprMax command and apply appropriate MPI configuration.

        Args:
            command: The command string to process.
            use_gprmax_mpi: Whether to use GprMax's built-in MPI.
            mpi_config: The MPIConfig object for mpirun wrapper mode.
            n_vcpus: Number of vCPUs requested by user.
            cmd_idx: Index of the command in the commands list (for error msgs).

        Returns:
            Either a string command or a Command object with MPI config.

        Raises:
            ValueError: If the command has incorrect MPI configuration.
        """
        has_mpi_flag = "-mpi" in command
        has_mpi_no_spawn = "--mpi-no-spawn" in command
        has_mpirun = "mpirun" in command

        if use_gprmax_mpi:
            # Mode 1: Using GprMax's built-in MPI
            if has_mpirun:
                raise ValueError(
                    f"Command {cmd_idx}: Cannot use 'mpirun' when "
                    f"use_gprmax_mpi=True. GprMax's built-in MPI uses the "
                    f"-mpi flag instead.\n"
                    f"Command: {command}\n"
                    f"Either:\n"
                    f"  1. Remove 'mpirun' and use the -mpi flag, or\n"
                    f"  2. Set use_gprmax_mpi=False to use mpirun wrapper mode")

            if has_mpi_no_spawn:
                raise ValueError(
                    f"Command {cmd_idx}: The --mpi-no-spawn flag is only for "
                    f"mpirun wrapper mode.\n"
                    f"Command: {command}\n"
                    f"When use_gprmax_mpi=True, use the -mpi flag instead.\n"
                    f"Example: 'python -m gprMax input.in -n 60 -mpi 61'")

            if not has_mpi_flag:
                raise ValueError(
                    f"Command {cmd_idx}: Missing -mpi flag when "
                    f"use_gprmax_mpi=True.\n"
                    f"Command: {command}\n"
                    f"GprMax's built-in MPI requires the -mpi flag.\n"
                    f"Example: 'python -m gprMax input.in -n 60 -mpi 61'\n"
                    f"Or set use_gprmax_mpi=False if you want to use mpirun "
                    f"wrapper mode.")

            # GprMax built-in MPI mode - return command as is
            return command

        elif mpi_config is not None:
            # Mode 2: Using mpirun wrapper (n_vcpus was specified)
            if has_mpi_flag and not has_mpi_no_spawn:
                raise ValueError(
                    f"Command {cmd_idx}: Found -mpi flag but --mpi-no-spawn "
                    f"is missing.\n"
                    f"Command: {command}\n"
                    f"When using mpirun wrapper (use_gprmax_mpi=False with "
                    f"n_vcpus specified), you must use --mpi-no-spawn flag.\n"
                    f"Example: 'python -m gprMax input.in -n 3 --mpi-no-spawn'\n"
                    f"Or set use_gprmax_mpi=True to use GprMax's built-in MPI "
                    f"with -mpi flag.")

            if has_mpirun:
                raise ValueError(
                    f"Command {cmd_idx}: Command already contains 'mpirun'.\n"
                    f"Command: {command}\n"
                    f"When n_vcpus is specified, the mpirun wrapper is added "
                    f"automatically. Please remove 'mpirun' from your command.\n"
                    f"Example: 'python -m gprMax input.in -n 3 --mpi-no-spawn'")

            if not has_mpi_no_spawn:
                raise ValueError(
                    f"Command {cmd_idx}: Missing --mpi-no-spawn flag.\n"
                    f"Command: {command}\n"
                    f"When using mpirun wrapper mode (use_gprmax_mpi=False "
                    f"with n_vcpus={n_vcpus}), the --mpi-no-spawn flag is "
                    f"required.\n"
                    f"Example: 'python -m gprMax input.in -n 3 --mpi-no-spawn'\n"
                    f"See: https://docs.gprmax.com/en/latest/openmp_mpi.html")

            # mpirun wrapper mode - wrap with Command and MPI config
            return Command(command, mpi_config=mpi_config)

        else:
            # Mode 3: No MPI (n_vcpus not specified)
            if has_mpi_flag or has_mpi_no_spawn or has_mpirun:
                raise ValueError(
                    f"Command {cmd_idx}: MPI-related flags found but MPI is "
                    f"not configured.\n"
                    f"Command: {command}\n"
                    f"The command contains MPI flags (-mpi, --mpi-no-spawn, "
                    f"or mpirun) but n_vcpus was not specified.\n"
                    f"Either:\n"
                    f"  1. Specify n_vcpus to enable MPI, or\n"
                    f"  2. Remove MPI flags if you want to run without MPI\n"
                    f"Example without MPI: "
                    f"'python -m gprMax input.in -n 60'")

            # No MPI mode - return command as is
            return command
