"""SCHISM module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


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
        self.simulator = "schism"

    def run(self,
            input_dir: str,
            *,
            on: types.ComputationalResources,
            num_scribes: int = 1,
            storage_dir: Optional[str] = "",
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            resubmit_on_preemption: bool = False,
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
        """
        return super().run(input_dir,
                           on=on,
                           num_scribes=num_scribes,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           use_hwthread=use_hwthread,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
