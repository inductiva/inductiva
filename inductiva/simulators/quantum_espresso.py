"""Class to run commands on Quantum Espresso."""
from typing import List, Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class QuantumEspresso(simulators.Simulator):
    """Class to run commands on Quantum Espresso."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the Quantum Espresso simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "quantumespresso"

    @property
    def name(self):
        """Get the name of the this simulator."""
        return "Quantum-Espresso"

    def run(self,
            input_dir: str,
            commands: List[str],
            *,
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            storage_dir: Optional[str] = "",
            on: types.ComputationalResources,
            extra_metadata: Optional[dict] = None,
            resubmit_on_preemption: bool = False,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            n_vcpus: Number of vCPUs to use in the mpi commands. If not provided
                (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available.
            on: The computatißonal resource to launch the simulation on.
            storage_dir: Parent directory for storing simulation
                               results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
        """
        return super().run(input_dir,
                           on=on,
                           n_vcpus=n_vcpus,
                           commands=commands,
                           storage_dir=storage_dir,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           container_image=self._image_uri,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
