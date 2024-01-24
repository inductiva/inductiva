"""OpenFOAM module of the API for fluid dynamics."""
from typing import Optional

from inductiva import types, tasks, simulators

AVAILABLE_OPENFOAM_VERSIONS = ["foundation", "esi"]


@simulators.simulator.mpi_enabled
class OpenFOAM(simulators.Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API.

    Users can choose between the ESI or the Foundation version
    by selecting the version on the initiliasation. Be aware, that
    some input files may only work for a specific version.
    """

    def __init__(self, version: str = "foundation"):
        if version not in AVAILABLE_OPENFOAM_VERSIONS:
            raise ValueError("Version not currently supported."
                             f"Available: {AVAILABLE_OPENFOAM_VERSIONS}")

        super().__init__()
        self.api_method_name = f"fvm.openfoam_{version}.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: types.Commands,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        extra_metadata: Optional[dict] = None,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            commands: List of commands to run using the OpenFOAM simulator.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata)
