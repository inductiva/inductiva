"""FEniCSx module of the API for finite element analysiss."""

from typing import Optional, List
from uuid import UUID

from inductiva import types, tasks
from inductiva.simulation import Simulator


class FEniCSx(Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self, api_method: str = "fem.fenicsx.run_simulation"):
        super().__init__()
        self.api_method = api_method

    @property
    def api_method_name(self) -> str:
        return self.api_method

    def run(
        self,
        mesh_path: types.Path,
        bc_path: types.Path,
        material_path: types.Path,
        resource_pool_id: Optional[UUID] = None,
        n_cores: int = 1,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            mesh_path: Path to the mesh file.
            bc_path: Path to the boundary conditions file.
            material_path: Path to the material file.
            n_cores: Number of MPI cores to use for the simulation.
            other arguments: See the documentation of the base class.
        """
        return super().run(mesh_path,
                           bc_path,
                           material_path,
                           resource_pool_id=resource_pool_id,
                           n_cores=n_cores,
                           run_async=run_async)
