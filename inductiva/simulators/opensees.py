"""OpenSees module of the API."""
from typing import List, Optional
import logging

from inductiva import types, tasks, simulators
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig

AVAILABLE_OPENSEES_INTERFACES = ["python", "tcl"]


@simulators.simulator.mpi_enabled
class OpenSees(simulators.Simulator):
    """Class to invoke a generic OpenSees simulation on the API."""

    def __init__(self,
                 /,
                 version: Optional[str] = None,
                 use_dev: bool = False,
                 interface: str = "python"):
        """Initialize the OpenSees simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
            interface (str): The interface to use for interacting with the simulator.
                Can be either "python" (default) or "tcl". 
        """
        if interface.lower() not in AVAILABLE_OPENSEES_INTERFACES:
            raise ValueError(
                f"Interface '{interface}' for OpenSees is not supported. "
                f"Available interfaces are: "
                f"{AVAILABLE_OPENSEES_INTERFACES}")
        self._interface = interface.lower()
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "opensees"

    @property
    def name(self):
        """Get the name of the simulator."""
        if self._interface == "python":
            return "OpenSeesPy"
        else:
            return "OpenSees"

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
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
            on: The computational resource to launch the simulation on.
            sim_config_filename: Name of the simulation configuration file.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
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
        
        if self._version == "2.5.0":
            logging.warning(
                "Opensees version 2.5.0 does not support parallel execution."
                "Changing to `n_vcpus` to 1.")
            n_vcpus=1
            if self._interface == "python":
                logging.warning(
                "Opensees version 2.5.0 does not support `python` as"
                " an interface. Changing to `tcl`.")

        if n_vcpus > on.

        if self._interface == "python":

            commands = [f"python {sim_config_filename}"]

            #If we use any mpi flag we will use mpi (otherwise run python only)
            if n_vcpus is not None or use_hwthread is True:
                logging.info("\nMPI flag detected (n_vcpus or use_hwthread).\n"
                             "We are going to run your python file with mpirun "
                             "(Parallel Execution).")
                mpi_config = MPIConfig(version="4.1.6",
                                       use_hwthread_cpus=use_hwthread,
                                       **({
                                           "np": n_vcpus
                                       } if n_vcpus is not None else {}))

                commands = [
                    Command(f"python {sim_config_filename}",
                            mpi_config=mpi_config)
                ]
            else:
                logging.info(
                    "\nNo MPI flag detected (n_vcpus or use_hwthread).\n"
                    "We are going to run your python file with python "
                    "(Sequential Execution).")
        else:
            mpi_config = MPIConfig(version="4.1.6",
                                   use_hwthread_cpus=use_hwthread,
                                   **({
                                       "np": n_vcpus
                                   } if n_vcpus is not None else {}))

            commands = [
                Command(f"OpenSeesMP {sim_config_filename}",
                        mpi_config=mpi_config)
            ]
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
