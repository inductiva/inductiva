"""OpenFOAM module of the API for fluid dynamics."""
from typing import Optional

from inductiva import types, tasks, simulators

AVAILABLE_OPENFOAM_BRANCHES = ["foundation", "esi"]


@simulators.simulator.mpi_enabled
class OpenFOAM(simulators.Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API.

    Users can choose between the ESI or the Foundation version
    by selecting the version on the initiliasation. Be aware, that
    some input files may only work for a specific version.
    """

    def __init__(self,
                 /,
                 branch: str = "foundation",
                 version: Optional[str] = None,
                 use_dev: bool = False):
        """Initialize the OpenFOAM simulator.

        Args:
            branch (str): The branch of OpenFOAM to use. Available branches are
                "foundation" and "esi". Default is "foundation".
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        if branch not in AVAILABLE_OPENFOAM_BRANCHES:
            raise ValueError(
                f"Branch {branch} of OpenFOAM is not supported. "
                f"Available branches are: {AVAILABLE_OPENFOAM_BRANCHES}")

        self._branch = branch
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = f"fvm.openfoam_{branch}.run_simulation"

    def _get_simulator_name(self):
        """Get the name of the simulator."""
        return "OpenFOAM-" + self._branch

    def run(self,
            input_dir: types.Path,
            commands: types.Commands,
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            extra_metadata: Optional[dict] = None,
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
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           n_vcpus=n_vcpus,
                           use_hwthread=use_hwthread,
                           extra_metadata=extra_metadata,
                           **kwargs)
