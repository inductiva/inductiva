"""Octopus module of the API."""
from typing import Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class Octopus(simulators.Simulator):
    """Class to invoke a generic Octopus simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the Octopus simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
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
            other arguments: See the documentation of the base class.
        """

        on.

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           **kwargs)
