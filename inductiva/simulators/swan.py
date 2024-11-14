"""SWAN module of the API."""
from typing import List, Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class SWAN(simulators.Simulator):
    """Class to invoke a generic SWAN simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SWAN simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "swan"

    def run(
        self,
        input_dir: Optional[str],
        sim_config_filename: Optional[str] = None,
        *,
        remote_assets: Optional[List[str]] = None,
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
        """

        if command not in ("swanrun", "swan.exe"):
            raise ValueError("Invalid command. Use 'swanrun' or 'swan.exe'.")

        if sim_config_filename is None and command == "swanrun":
            raise ValueError("Simulation configuration file "
                             "(sim_config_filename) not provided.\n"
                             "When using 'swanrun' it is mandatory to provide "
                             "sim_config_filename.")

        return super().run(input_dir,
                           on=on,
                           input_filename=sim_config_filename,
                           storage_dir=storage_dir,
                           command=command,
                           n_vcpus=n_vcpus,
                           use_hwthread=use_hwthread,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
