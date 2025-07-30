"""Octopus module of the API."""
import logging
from typing import Literal, Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class Octopus(simulators.Simulator):
    """Class to invoke a generic Octopus simulation on the API."""

    def __init__(self,
                 /,
                 version: Optional[str] = None,
                 use_dev: bool = False,
                 device: Literal["auto", "cpu", "gpu"] = "auto"):
        """Initialize the Octopus simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
            device (str): The device to use for the simulation. If "auto",
                the device will be selected automatically. If "cpu", the
                simulation will run on the CPU. If "gpu", the simulation
                will run on the GPU.
        """
        super().__init__(version=version, use_dev=use_dev, device=device)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "octopus"

    def run(self,
            input_dir: Optional[str],
            *,
            commands: types.Commands,
            on: types.ComputationalResources,
            n_vcpus: Optional[int] = None,
            n_mpi_processes: Optional[int] = None,
            n_omp_threads: Optional[int] = None,
            use_hwthread: bool = True,
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
            commands: List of commands to run.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
            n_omp_threads: Number of OpenMP threads to use in the simulation.
                If not provided, it defaults to the number of available vcpus
                of the machine group where the task will run.
            n_mpi_processes: Number of MPI processes that will run the
                simulation. If not provided, it defaults to 1.
            use_hwthread: If specified, Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available.
            storage_dir: Directory for storing simulation results.
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

        available_mpi_slots = on.get_available_mpi_slots(
            use_hwthread=use_hwthread)

        available_vcpus = on.available_vcpus

        if n_vcpus is None:
            logging.info(
                "Param n_vcpus not set. Defaulting to the number of "
                "available vcpus (%s).\n", available_vcpus)
            n_vcpus = available_vcpus

        if n_mpi_processes is None:
            n_mpi_processes = 1

            logging.info("Param n_mpi_processes not set. Defaulting to %s.\n",
                         n_mpi_processes)

        if n_omp_threads is None:
            n_omp_threads = max(n_vcpus // n_mpi_processes, 1)
            logging.info("Param n_omp_threads not set. Defaulting to %s.\n",
                         n_omp_threads)

        if n_mpi_processes > available_mpi_slots:
            raise ValueError(
                f"n_mpi_processes ({n_mpi_processes}) exceeds the number of "
                f"MPI slots available on the machine ({available_mpi_slots})")

        if n_vcpus <= 0:
            raise ValueError("n_vcpus must be larger than 0")
        if n_mpi_processes <= 0:
            raise ValueError("n_mpi_processes must be larger than 0")
        if n_omp_threads <= 0:
            raise ValueError("n_omp_threads must be larger than 0")

        requested_vcpus = n_omp_threads * n_mpi_processes
        if requested_vcpus > available_vcpus:
            raise ValueError(
                f"n_mpi_processes * n_omp_threads ({n_mpi_processes} * "
                f"{n_omp_threads} = {requested_vcpus}) can't be larger than "
                f"{available_vcpus} (the number of available VCPUs in "
                "the specified machine).\n"
                "Update `n_mpi_processes` or `n_omp_threads`.")

        mpi_config = None
        if n_mpi_processes > 1:
            mpi_kwargs = {}
            mpi_kwargs["use_hwthread_cpus"] = use_hwthread
            mpi_kwargs["np"] = n_mpi_processes
            mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        # Will add the mpiconfig to the octopus command if needed
        processed_commands = []
        for command in commands:
            is_string = isinstance(command, str)
            # The only command to run in parallel
            is_octopus = "octopus" in command
            already_mpi = "mpirun" in command

            if (is_string and is_octopus and not already_mpi):
                processed_commands.append(
                    Command(
                        "octopus ",
                        mpi_config=mpi_config,
                        env={"OMP_NUM_THREADS": str(n_omp_threads)},
                    ))
            else:
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
