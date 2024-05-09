"""AmrWind module of the API for numerical simulations of fluid flows."""
from typing import Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class AmrWind(simulators.Simulator):
    """Class to invoke a generic AmrWind simulation on the API.
    """

    def __init__(self):

        super().__init__()
        self.api_method_name = "amrWind.amrWind.run_simulation"

    def run(self,
            input_dir: types.Path,
            sim_config_filename: str,
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            extra_metadata: Optional[dict] = None,
            storage_dir: Optional[types.Path] = "",
            on: Optional[types.ComputationalResources] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory of the simulation input files.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           on=on,
                           n_vcpus=n_vcpus,
                           storage_dir=storage_dir,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           input_filename=sim_config_filename,
                           **kwargs)
