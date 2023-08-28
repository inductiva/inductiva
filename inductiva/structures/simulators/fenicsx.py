"""FEniCSx module of the API for Finite Element Analysis."""

from typing import Optional
from uuid import UUID

import inductiva


class FEniCSx(inductiva.simulation.Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self, api_method: str = "fem.fenicsx.run_simulation"):
        super().__init__()
        self.api_method = api_method

    @property
    def api_method_name(self) -> str:
        return self.api_method

    def run(
        self,
        input_dir: inductiva.types.Path,
        mesh_filename: str,
        bcs_filename: str,
        material_filename: str,
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
    ) -> inductiva.tasks.Task:
        """Run the simulation.

        Args:
            mesh_filename: Mesh filename.
            bcs_filename: Boundary conditions filename.
            material_filename: Material filename.
            resource_pool_id: Optional UUID of the resource pool to use.
            run_async: Whether to run the simulation asynchronously.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           resource_pool_id=resource_pool_id,
                           run_async=run_async,
                           mesh_filename=mesh_filename,
                           bcs_filename=bcs_filename,
                           material_filename=material_filename)
