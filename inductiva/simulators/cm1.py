"""CM1 module of the API."""
from typing import List, Optional

from inductiva import types, tasks, simulators
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class CM1(simulators.Simulator):
    """Class to invoke a generic CM1 simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the CM1 simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "cm1"

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
            mode: str = "mpi",
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            sim_config_filename: Optional[str] = None,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            mode (str): The execution mode, either 'mpi' or 'openmp'.
            on: The computational resource to launch the simulation on.
            sim_config_filename: Name of the simulation configuration file.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
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
            other arguments: See the documentation of the base class.
        """
        mode = mode.lower()
        executable = "cm1.exe" if mode == "mpi" else "cm1_openmp.exe"
        mpi_config = None

        if mode == "mpi":
            mpi_kwargs = {"use_hwthread_cpus": use_hwthread}
            if n_vcpus is not None:
                mpi_kwargs["np"] = n_vcpus
            mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        commands = [
            Command(f"{executable} {sim_config_filename}",
                    mpi_config=mpi_config)
        ]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
