"""OpenFAST module of the API for wind turbine simulations."""
from typing import Optional

from inductiva import types, tasks, simulators


class OpenFAST(simulators.Simulator):
    """Class to invoke a generic OpenFAST simulation on the API.

    """

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the OpenFAST simulator.
        
        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = "openfast.openfast.run_simulation"

    def run(self,
            input_dir: str,
            commands: types.Commands,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            resubmit_on_preemption: bool = False,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            commands: List of commands to run using the OpenFAST simulator.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            other arguments: See the documentation of the base class.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
