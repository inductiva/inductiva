"""GROMACS module of the API"""

import pathlib
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
        output_dir: Optional[types.Path] = None,
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
    ) -> pathlib.Path:
        """Run a list of GROMACS commands.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run using the GROMACS simulator.
            output_dir: Path to the directory where the output files will be
              stored when running synchronously.
            resource_pool_id: UUID of the resource pool to use for the
              simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        return super().run(input_dir,
                           output_dir=output_dir,
                           resource_pool_id=resource_pool_id,
                           commands=commands,
                           run_async=run_async)
