"""GROMACS module of the API"""

from typing import Optional, List
from uuid import UUID

from inductiva import types
from inductiva.simulation import Simulator


class GROMACS(Simulator):
    """Class to invoke any GROMACS command on the API."""

    @property
    def api_method_name(self) -> str:
        return "md.gromacs.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: List[dict],
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
    ):
        """Run a list of GROMACS commands.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run using the GROMACS simulator.
            resource_pool_id: UUID of the resource pool to use for the
              simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        return super().run(input_dir,
                           resource_pool_id=resource_pool_id,
                           commands=commands,
                           run_async=run_async)
