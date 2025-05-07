""" SNL SWAN module of the API."""
from typing import List, Optional
from pathlib import Path

from inductiva import types, tasks, simulators
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class SNLSWAN(simulators.Simulator):
    """Class to invoke a generic SNL SWAN simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SNL SWAN simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "snlswan"

    @property
    def name(self):
        """Get the name of the this simulator."""
        return "Snl-Swan"

    def run(
        self,
        input_dir: Optional[str],
        *,
        sim_config_filename: Optional[str],
        remote_assets: Optional[List[str]] = None,
        project: Optional[str] = None,
        resubmit_on_preemption: bool = False,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        n_vcpus: Optional[int] = None,
        use_hwthread: bool = True,
        command: str = "swanrun",
        **kwargs,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
                Mandatory when using 'swanrun' command.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            on: The computational resource to launch the simulation on.
            storage_dir: Directory for storing simulation results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            command: The command to run the simulation. Default is 'swanrun'.
                The user can also specify 'swan.exe'.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
        """

        if command not in ("swanrun", "swan.exe"):
            raise ValueError("Invalid command. Use 'swanrun' or 'swan.exe'.")

        if sim_config_filename is None and command == "swanrun":
            raise ValueError("Simulation configuration file "
                             "(sim_config_filename) not provided.\n"
                             "When using 'swanrun' it is mandatory to provide "
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

        # Swanrun uses internal MPI
        # we call apptainer run ... swanrun ... -mpi np
        if command == "swanrun":

            machinefile_command = Command(
                "dd if=/dev/stdin of=machinefile",
                f"localhost slots={on.available_vcpus}")

            commands.append(machinefile_command)

            #if the user does not provide n_vcpus use all available by default
            mpi_flag = f"-mpi {n_vcpus or on.available_vcpus}"

            swanrun_command = Command(
                f"swanrun -input {config_file_only} {mpi_flag}")
            commands.append(swanrun_command)

        # we call mpirun ... apptainer ... Swan.exe
        # works with clusters
        elif command == "swan.exe":

            mpi_kwargs = {}
            #If the user does not provide n_vcpus mpi will use all available
            if n_vcpus is not None:
                mpi_kwargs["np"] = n_vcpus
            mpi_kwargs["use_hwthread_cpus"] = use_hwthread

            mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)
            swan_exe_command = Command(f"swan.exe {sim_config_filename}",
                                       mpi_config=mpi_config)
            commands.append(swan_exe_command)
        return super().run(input_dir,
                           on=on,
                           storage_dir=storage_dir,
                           commands=commands,
                           run_subprocess_dir=str(working_dir),
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
