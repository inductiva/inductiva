"""SWASH module of the API."""
from pathlib import Path
from typing import List, Optional

from inductiva import types, tasks, simulators
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class SWASH(simulators.Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SWASH simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "swash"

    def run(self,
            input_dir: Optional[str],
            sim_config_filename: str,
            *,
            command: str = "swashrun",
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            storage_dir: Optional[str] = "",
            on: types.ComputationalResources,
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            project: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            command: The command to run the simulation. Default is 'swashrun'.
                The user can also specify 'swash.exe'.
            storage_dir: Directory for storing simulation results.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
        """

        if command not in ("swashrun", "swash.exe"):
            raise ValueError("Invalid command. Use 'swashrun' or 'swash.exe'.")

        if sim_config_filename is None and command == "swashrun":
            raise ValueError("Simulation configuration file "
                             "(sim_config_filename) not provided.\n"
                             "When using 'swashrun' it is mandatory to provide "
                             "sim_config_filename.")

        commands = []

        path_config_filename = Path(sim_config_filename)

        if path_config_filename.is_absolute():
            raise ValueError("sim_config_filename must be a path relative to "
                             "the input directory.")

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename=sim_config_filename)

        working_dir = path_config_filename.parent
        config_file_only = path_config_filename.name

        # swashrun uses internal MPI
        # we call apptainer run ... swashrun ... -mpi np
        if command == "swashrun":

            machinefile_command = Command(
                "dd if=/dev/stdin of=machinefile",
                f"localhost slots={on.available_vcpus}")

            commands.append(machinefile_command)

            #if the user does not provide n_vcpus use all available by default
            mpi_flag = f"-mpi {n_vcpus or on.available_vcpus}"

            swashrun_command = Command(
                f"swashrun -input {config_file_only} {mpi_flag}")
            commands.append(swashrun_command)

        # we call mpirun ... apptainer ... swash.exe
        # works with clusters
        elif command == "swash.exe":

            mpi_kwargs = {}
            #If the user does not provide n_vcpus mpi will use all available
            if n_vcpus is not None:
                mpi_kwargs["np"] = n_vcpus
            mpi_kwargs["use_hwthread_cpus"] = use_hwthread

            mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)
            swash_exe_command = Command(f"swash.exe {sim_config_filename}",
                                        mpi_config=mpi_config)
            commands.append(swash_exe_command)

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           run_subprocess_dir=str(working_dir),
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
