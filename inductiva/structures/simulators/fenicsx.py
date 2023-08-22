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
        mesh_filename: str,
        bcs_filename: str,
        material_filename: str,
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
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
