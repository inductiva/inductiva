"""OpenFOAM module of the API for fluid dynamics."""
from typing import Optional

from inductiva import types, tasks, simulators

AVAILABLE_OPENFOAM_DISTRIBUTIONS = ["foundation", "esi"]


@simulators.simulator.mpi_enabled
class OpenFOAM(simulators.Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API.

    Users can choose between the ESI or the Foundation version
    by selecting the version on the initiliasation. Be aware, that
    some input files may only work for a specific version.
    """

    def __init__(self,
                 /,
                 distribution: str = "foundation",
                 version: Optional[str] = None,
                 use_dev: bool = False):
        """Initialize the OpenFOAM simulator.

        Args:
            distribution (str): The distribution of OpenFOAM to use. Available
                distributions are "foundation" and "esi". Default is
                "foundation".
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        if distribution not in AVAILABLE_OPENFOAM_DISTRIBUTIONS:
            raise ValueError(
                f"Distribution '{distribution}' of OpenFOAM is not supported. "
                f"Available distributions are: "
                f"{AVAILABLE_OPENFOAM_DISTRIBUTIONS}")

        self._distribution = distribution
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = f"fvm.openfoam_{distribution}.run_simulation"

    @property
    def name(self):
        """Get the name of the simulator."""
        return "OpenFOAM-" + self._distribution

    def run(self,
            input_dir: str,
            commands: types.Commands,
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            resubmit_on_preemption: bool = False,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            commands: List of commands to run using the OpenFOAM simulator.
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
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
