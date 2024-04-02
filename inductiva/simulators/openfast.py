"""Openfast module of the API for wind turbine simulations."""
from typing import Optional

from inductiva import types, tasks, simulators


class Openfast(simulators.Simulator):
    """Class to invoke a generic Openfast simulation on the API.

    """

    def __init__(self):

        super().__init__()
        self.api_method_name = "openfast.openfast.run_simulation"

    def run(self,
            input_dir: types.Path,
            commands: types.Commands,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            extra_metadata: Optional[dict] = None,
            provider_id: str = "GCP",
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            commands: List of commands to run using the Openfast simulator.
            number of hardware threads on the node, and use that as the
            number of slots available.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           n_vcpus=1,
                           use_hwthread=False,
                           provider_id=provider_id,
                           extra_metadata=extra_metadata)
