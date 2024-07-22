"""SWASH module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class SWASH(simulators.Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SWASH simulator.
        
        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = "sw.swash.run_simulation"

    def run(self,
            input_dir: str,
            sim_config_filename: str,
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            resubmit_on_preemption: bool = False,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
            storage_dir: Directory for storing simulation results.
        """
        return super().run(input_dir,
                           on=on,
                           input_filename=sim_config_filename,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
