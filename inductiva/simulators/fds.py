"""FDS simulator module of the API."""

from typing import Optional, Union
import logging

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


class FDS(simulators.Simulator):
    """Class to invoke a generic FDS simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the FDS simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "fds"

    def run(self,
            input_dir: Optional[str],
            sim_config_filename: str,
            *,
            on: types.ComputationalResources,
            n_vcpus: Optional[int] = None,
            n_omp_threads: Optional[int] = None,
            n_mpi_processes: Optional[int] = None,
            use_hwthread: bool = True,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            sim_config_filename: Name of the simulation configuration file.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
            n_omp_threads: Number of OpenMP threads to use in the simulation.
                If not provided, it defaults to the number of available vcpus
                of the machine group where the task will run.
            n_mpi_processes: Number of MPI processes that will run the
                simulation. If not provided, it defaults to 1. Note that the
                number of MPI processes can't exceed the number of meshes
                of the simulation case.
            use_hwthread: If specified Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            other arguments: See the documentation of the base class.
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
        """

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename=sim_config_filename)

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

            logging.info(
                "Param n_mpi_processes not set. Defaulting to %s. "
                "Note that the number of MPI processes for FDS "
                "simulations is limited by the number of meshes "
                "in the simulation case.\n", n_mpi_processes)

        if n_omp_threads is None:
            n_omp_threads = max(n_vcpus // n_mpi_processes, 1)
            logging.info("Param n_omp_threads not set. Defaulting to %s.\n",
                         n_omp_threads)

        if n_mpi_processes > available_mpi_slots:
            raise ValueError(
                f"n_mpi_processes ({n_mpi_processes}) exceeds the number of "
                f"MPI slots available on the machine ({available_mpi_slots})")

        if n_vcpus == 0:
            raise ValueError("n_vcpus must be larger than 0")
        if n_mpi_processes == 0:
            raise ValueError("n_mpi_processes must be larger than 0")
        if n_omp_threads == 0:
            raise ValueError("n_omp_threads must be larger than 0")

        requested_vcpus = n_omp_threads * n_mpi_processes
        if requested_vcpus > available_vcpus:
            raise ValueError(
                f"n_mpi_processes * n_omp_threads ({n_mpi_processes} * "
                f"{n_omp_threads} = {requested_vcpus}) can't be larger than "
                f"{available_vcpus} (the number of available VCPUs in "
                "the specified machine)")

        mpi_config = None
        if n_mpi_processes > 1:
            mpi_kwargs = {}
            mpi_kwargs["use_hwthread_cpus"] = use_hwthread
            mpi_kwargs["np"] = n_mpi_processes
            mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        commands = [
            Command(
                "/opt/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux "
                f"{sim_config_filename}",
                mpi_config=mpi_config,
                env={"OMP_NUM_THREADS": str(n_omp_threads)},
            )
        ]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           **kwargs)
