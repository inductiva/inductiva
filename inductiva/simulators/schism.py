"""SCHISM module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class SCHISM(simulators.Simulator):
    """Class to invoke a generic SCHISM simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "schism.schism.run_simulation"

    def run(self,
            input_dir: types.Path,
            num_scribes: int = 1,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            use_hwthread: bool = True,
            extra_metadata: Optional[dict] = None,
            n_vcpus: int = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory of the simulation input files.
            num_scribes: The num_scribes as per the simulator documentation.
            # pylint: disable=line-too-long
            https://schism-dev.github.io/schism/master/getting-started/running-model.html
            # pylint: enable=line-too-long
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Directory for storing simulation results.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            n_vcpus: Number of virtual cpus
        """
        return super().run(input_dir,
                           on=on,
                           num_scribes=num_scribes,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           **kwargs)
