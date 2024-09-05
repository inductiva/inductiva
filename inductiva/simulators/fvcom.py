"""FVCOM simulator module of the API."""
from typing import Optional

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class FVCOM(simulators.Simulator):
    """Class to invoke a generic FVCOM simulation on the API.

    """

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the FVCOM simulator.
        
        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = "fvcom.fvcom.run_simulation"

    def run(self,
            input_dir: str,
            case_name: str,
            debug: int = 0,
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            run_from: Optional[str] = None,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            extra_metadata: Optional[dict] = None,
            on: Optional[types.ComputationalResources] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            casename: Name of the simulation case.
            debug: Debug level of the simulation (from 0 to 7).
            run_from: Path (relative to the input directory) to the directory
                where the simulation nml file is located. If not provided, the
                input directory is used.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
            other arguments: See the documentation of the base class.
        """
        if debug < 0 or debug > 7:
            raise ValueError("Debug level must be between 0 and 7.")
        
        command = f"./fvcom --CASENAME={case_name} --dbg={debug}"
        return super().run(input_dir,
                           on=on,
                           command=command,
                           n_vcpus=n_vcpus,
                           run_from=run_from,
                           storage_dir=storage_dir,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
