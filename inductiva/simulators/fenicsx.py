"""FEniCSx module of the API for Finite Element Analysis."""

from typing import List, Optional

from inductiva import simulators, types, tasks


class FEniCSx(simulators.Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the FEniCSx simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "fenicsx"

    def _get_image_uri(self):
        return None

    def run(self,
            input_dir: Optional[str],
            geometry_filename: str,
            bcs_filename: str,
            material_filename: str,
            *,
            on: types.ComputationalResources,
            mesh_filename: Optional[str] = None,
            mesh_info_filename: Optional[str] = None,
            storage_dir: Optional[str] = "",
            remote_assets: Optional[List[str]] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            on: The computational resource to launch the simulation on.
            geometry_filename: Geometry filename.
            bcs_filename: Boundary conditions filename.
            material_filename: Material filename.
            mesh_filename: Mesh filename.
            mesh_info_filename: Mesh information filename.
            storage_dir: Parent directory for storing simulation results.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
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
                           mesh_filename=mesh_filename,
                           bcs_filename=bcs_filename,
                           storage_dir=storage_dir,
                           on=on,
                           remote_assets=remote_assets,
                           **kwargs)
