"""OpenSees module of the API."""
import logging
from typing import Literal, Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig

AVAILABLE_OPENSEES_INTERFACES = ["python", "tcl", "eesd"]


@simulators.simulator.mpi_enabled
class OpenSees(simulators.Simulator):
    """Class to invoke a generic OpenSees simulation on the API."""

    def __init__(self,
                 /,
                 version: Optional[str] = None,
                 use_dev: bool = False,
                 interface: Literal["python", "tcl", "eesd"] = "python"):
        """Initialize the OpenSees simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
            interface (str): The interface to use for interacting with the
                simulator. Can be either "python" (default) or "tcl".
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
        elif self._interface == "eesd":
            return "OpenSees-EESD"
        else:
            return "OpenSees"

    def _build_commands(self, sim_config_filename, mpi_config):
        """ Creates the list of commands based on the sim_config_filename and
        the provided mpi_config.
        """
        if self._interface == "python":
            return [
                Command(f"python {sim_config_filename}", mpi_config=mpi_config)
            ]

        if self._interface == "eesd":
            return [Command(f"OpenSees {sim_config_filename}")]

        return [
            Command(f"OpenSeesMP {sim_config_filename}", mpi_config=mpi_config)
        ]

    def _build_mpi_config(self, n_vcpus: int, use_hwthread: bool) -> MPIConfig:
        """ Creates an MPIConfig based on n_vcpus and use_hwthread.
        """
        return MPIConfig(version="4.1.6",
                         use_hwthread_cpus=use_hwthread,
                         **({
                             "np": n_vcpus
                         } if n_vcpus is not None else {}))

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
            sim_config_filename: Optional[str],
            n_vcpus: Optional[int] = None,
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

        if self._version == "2.5.0":
            logging.warning(
                "Opensees version 2.5.0 does not support parallel execution."
                "Changing to `n_vcpus` to 1.")
            n_vcpus = 1
            if self._interface == "python":
                raise ValueError(
                    "Opensees version 2.5.0 does not support `python` as"
                    " an interface. Please use `interface='tcl'`.")

        self._check_vcpus(n_vcpus, on)

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename=sim_config_filename)

        mpi_config = self._build_mpi_config(n_vcpus, use_hwthread)

        if self._interface == "python":
            if n_vcpus or use_hwthread:
                logging.info("\nMPI flag detected (n_vcpus or use_hwthread)."
                             "\nThe simulation will run with mpirun "
                             "(Parallel Execution).")
                commands = self._build_commands(sim_config_filename, mpi_config)
            else:
                logging.info(
                    "\nNo MPI flag detected (n_vcpus or use_hwthread)."
                    "\nThe simulation will run sequentially with python.")
                commands = [f"python {sim_config_filename}"]
        else:
            commands = self._build_commands(sim_config_filename, mpi_config)

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           on_finish_cleanup=on_finish_cleanup,
                           **kwargs)
