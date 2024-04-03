"""OpenFAST module of the API for wind turbine simulations."""
from typing import Optional

from inductiva import types, tasks, simulators


class OpenFAST(simulators.Simulator):
    """Class to invoke a generic OpenFAST simulation on the API.

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
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            commands: List of commands to run using the OpenFAST simulator.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           **kwargs)
