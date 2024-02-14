"""SWASH module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators
from inductiva.utils import meta


@simulators.simulator.mpi_enabled
class SWASH(simulators.Simulator):
    """Class to invoke a generic SWASH simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sw.swash.run_simulation"

    @meta.deprecated_arg(n_cores="n_vcpus")
    def run(self,
            input_dir: types.Path,
            sim_config_filename: str,
            n_vcpus: Optional[int] = None,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            extra_metadata: Optional[dict] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            n_vcpus: Number of vCPUs (all by default) to use for the simulation.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Directory for storing simulation results.
        """
        return super().run(input_dir,
                           on=on,
                           input_filename=sim_config_filename,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           extra_metadata=extra_metadata)
