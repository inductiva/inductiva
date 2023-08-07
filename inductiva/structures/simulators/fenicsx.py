"""FEniCSx module of the API for finite element analysiss."""

from typing import Optional
from uuid import UUID

from inductiva import types, tasks
from inductiva import simulation


class FEniCSx(simulation.Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self, api_method: str = "fem.fenicsx.run_simulation"):
        super().__init__()
        self.api_method = api_method

    @property
    def api_method_name(self) -> str:
        return self.api_method

    def run(
        self,
        input_dir: types.Path,
        geometry_filename: types.Path,
        mesh_filename: types.Path,
        bc_filename: types.Path,
        material_filename: types.Path,
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            geometry_filename: Filename of the geometry file.
            mesh_filename: Filename of the mesh file.
            bc_filename: Filename of the boundary conditions file.
            material_filename: Filename of the material file.
            resource_pool_id: Optional UUID of the resource pool to use.
            run_async: Whether to run the simulation asynchronously.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           geometry_filename,
                           mesh_filename,
                           bc_filename,
                           material_filename,
                           resource_pool_id=resource_pool_id,
                           run_async=run_async)
