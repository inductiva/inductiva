"""Reef3D simulator module of the API."""

from typing import Optional

from inductiva import simulators, types, tasks
from inductiva.utils import meta


@simulators.simulator.mpi_enabled
class REEF3D(simulators.Simulator):
    """Class to invoke a generic FDS simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "reef3d.reef3d.run_simulation"

    @meta.deprecated_arg(n_cores="n_vcpus")
    def run(self,
            input_dir: types.Path,
            n_vcpus: Optional[int] = None,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            extra_metadata: Optional[dict] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           on=on,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           extra_metadata=extra_metadata)
