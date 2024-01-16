"""XBeach module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class XBeach(simulators.Simulator):
    """Class to invoke a generic XBeach simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sw.xbeach.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        sim_config_filename: Optional[str] = "params.txt",
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        extra_metadata: Optional[dict] = None,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Directory for storing simulation results.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           input_filename=sim_config_filename,
                           on=on,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata)
