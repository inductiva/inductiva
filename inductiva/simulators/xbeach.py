"""XBeach module of the API."""
import logging
from typing import Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class XBeach(simulators.Simulator):
    """Class to invoke a generic XBeach simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the XBeach simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "xbeach"

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            sim_config_filename: Optional[str] = None,
            export_vtk: bool = False,
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
                Deprecated: This parameter is no longer used and will be removed
                in a future version. Please use `params.txt` in the input
                directory instead.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available.
            export_vtk (bool): If True, after the XBeach simulation finishes 
                generate the VTK file via `xbeach_animator --export-vtk`.
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

        if sim_config_filename is not None:
            logging.warning(
                "Deprecated: This parameter is no longer used and will be "
                "removed in a future version. Please use `params.txt` in the "
                "input directory instead.")

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename="params.txt")

        mpi_kwargs = {}
        mpi_kwargs["use_hwthread_cpus"] = use_hwthread
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus

        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)
        commands = [Command("xbeach", mpi_config=mpi_config)]

        # Conditionally append the VTK-export step
        if export_vtk:
            commands.append(
                Command("python3 /usr/local/bin/xbeach_animator "
                        "--input-dir . --output-dir . --export-vtk"))

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           **kwargs)
