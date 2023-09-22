"""FEniCSx module of the API for Finite Element Analysis."""

from typing import Optional

from inductiva import types, tasks, resources
from inductiva.simulators import Simulator


class FEniCSx(Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "fem.fenicsx.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        mesh_filename: str,
        bcs_filename: str,
        material_filename: str,
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            mesh_filename: Mesh filename.
            bcs_filename: Boundary conditions filename.
            material_filename: Material filename.
            machine_group: The machine group to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir,
                           machine_group=machine_group,
                           run_async=run_async,
                           mesh_filename=mesh_filename,
                           bcs_filename=bcs_filename,
                           material_filename=material_filename)
