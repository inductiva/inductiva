"""FEniCSx module of the API for Finite Element Analysis."""

from typing import Optional

from inductiva import simulators, types, tasks


class FEniCSx(simulators.Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "fem.fenicsx.run_simulation"

    def run(self,
            input_dir: types.Path,
            geometry_filename: str,
            bcs_filename: str,
            material_filename: str,
            mesh_filename: Optional[str] = None,
            mesh_info_filename: Optional[str] = None,
            on: Optional[types.ComputationalResources] = None,
            storage_dir: Optional[types.Path] = "",
            extra_metadata: Optional[dict] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            geometry_filename: Geometry filename.
            bcs_filename: Boundary conditions filename.
            material_filename: Material filename.
            mesh_filename: Mesh filename.
            mesh_info_filename: Mesh information filename.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Parent directory for storing simulation results.
            **kwargs: Arbitrary keyword arguments, including:
                - global_refinement_meshing_factor (float): Factor for global
                  mesh refinement.
                - local_refinement_meshing_factor (float): Factor for local mesh
                  refinement.
                - smoothing_parameter (float): Mesh smoothing parameter.
                - mesh_element_family (str): Mesh element family.
                - mesh_element_order (int): Order of the mesh element.
                - mesh_quadrature_rule (str): Mesh quadrature rule.
                - mesh_quadrature_degree (int): Mesh quadrature degree.
            other arguments: See the documentation of the base class.
        """

        return super().run(input_dir,
                           mesh_info_filename=mesh_info_filename,
                           material_filename=material_filename,
                           geometry_filename=geometry_filename,
                           extra_metadata=extra_metadata,
                           mesh_filename=mesh_filename,
                           bcs_filename=bcs_filename,
                           storage_dir=storage_dir,
                           on=on,
                           **kwargs)
