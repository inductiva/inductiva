"""FEniCSx module of the API for Finite Element Analysis."""

from typing import Optional

from inductiva import simulators, types, tasks


class FEniCSx(simulators.Simulator):
    """Class to invoke a generic FEniCSx simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "fem.fenicsx.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        geometry_filename: str,
        bcs_filename: str,
        material_filename: str,
        global_refinement_meshing_factor: float = 1.0,
        local_refinement_meshing_factor: float = 1.0,
        smoothing_meshing_parameter: float = 10.0,
        mesh_element_family: str = "CG",
        mesh_element_order: int = 1,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        extra_metadata: Optional[dict] = None,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            geometry_filename: Geometry filename.
            bcs_filename: Boundary conditions filename.
            material_filename: Material filename.
            global_refinement_meshing_factor (float): The refinement factor for
              global refinement of the mesh. A higher value results in a finer
              mesh overall, increasing the number of elements in the entire
              mesh, and leading to a more detailed representation of the
              geometry. Use this factor when you want to globally refine the
              mesh uniformly, without specific local focus.
            local_refinement_meshing_factor (float): The refinement factor for
              local refinement of the mesh. This factor controls the local
              refinement level of the mesh and is typically used for refining
              specific regions or features of the mesh. A higher value for this
              factor indicates a finer mesh in the regions of interest,
              providing more detailed resolution around certain features. Use
              this factor when you want to focus on refining specific areas
              while keeping the rest of the mesh less refined.
            smoothing_parameter (float): The smoothing parameter for mesh
              generation. It controls the amount of mesh smoothing applied to
              the generated mesh. Adjust this parameter for improved mesh
              quality.
            mesh_element_family (str): The type of mesh element family.
            mesh_element_order (int): The (polynomial) order of the mesh
              element.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Parent directory for storing simulation results.
            other arguments: See the documentation of the base class.
        """

        return super().run(
            input_dir,
            on=on,
            geometry_filename=geometry_filename,
            bcs_filename=bcs_filename,
            material_filename=material_filename,
            global_refinement_meshing_factor=global_refinement_meshing_factor,
            local_refinement_meshing_factor=local_refinement_meshing_factor,
            smoothing_meshing_parameter=smoothing_meshing_parameter,
            mesh_element_family=mesh_element_family,
            mesh_element_order=mesh_element_order,
            storage_dir=storage_dir,
            extra_metadata=extra_metadata,
        )
