"""SCHISM module of the API."""
from typing import List, Optional

from inductiva import types, tasks, simulators
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class SCHISM(simulators.Simulator):
    """Class to invoke a generic SCHISM simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SCHISM simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "schism"

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
            num_scribes: int = 1,
            storage_dir: Optional[str] = "",
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            project: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            num_scribes: The num_scribes as per the simulator documentation.
            https://schism-dev.github.io/schism/master/getting-started/running-model.html
            storage_dir: Directory for storing simulation results.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            n_vcpus: Number of virtual cpus
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
        """
        mpi_kwargs = {}
        mpi_kwargs["use_hwthread_cpus"] = use_hwthread
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus

        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        commands = [
            "mkdir -p outputs",
            Command(f"/schism/build/bin/pschism {num_scribes}",
                    mpi_config=mpi_config)
        ]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
