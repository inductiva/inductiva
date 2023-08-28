"""Utils to create the mesh."""

from typing import List, Optional

import gmsh

from . import geometry_utils


def add_mesh_field_distance(mesh_field_id: int, curves_list: List[int],
                            num_points: int) -> None:
    """Add a distance-based mesh size field to control element size in Gmsh.

    - Use gmsh.model.mesh.field.add("Distance", mesh_field_id) to add a
    distance-based mesh size field, where mesh_field_id is a unique identifier
    for the field.
    - Set the curves or points where you want to apply the distance-based field
    using gmsh.model.mesh.field.setNumbers(field_id, "CurvesList", curves_list)
    or gmsh.model.mesh.field.setNumbers(field_id, "PointsList", points_list). 
    Adjust the points_list of the field (i.e., how many points to use along
    curves or points) with 
    gmsh.model.mesh.field.setNumber(field_id, "Sampling", num_points).

    Args:
        mesh_field_id (int): A unique identifier for the distance-based
          mesh size field.
        curves_list (List[int]): List of curve IDs where the
          distance-based field will be applied.
        num_points (int): Number of points to use along curves to adjust
          the field size.
    """

    # Add a distance-based mesh size field using the given field ID
    gmsh.model.mesh.field.add("Distance", mesh_field_id)

    # Set the curves where the distance-based field will be applied
    gmsh.model.mesh.field.setNumbers(mesh_field_id, "CurvesList", curves_list)

    # Adjust the points_list of the field to control mesh size along
    # curves
    gmsh.model.mesh.field.setNumber(mesh_field_id, "Sampling", num_points)


def add_mesh_field_threshold(mesh_field_id: int, ref_mesh_field_id: int,
                             size_min: float, size_max: float, dist_min: float,
                             dist_max: float) -> None:
    """Add a threshold-based mesh size field to control element size in Gmsh.

    - Use gmsh.model.mesh.field.add("Threshold", mesh_field_id) to add a 
    threshold-based mesh size field, where mesh_field_id is a unique identifier
    for the field.
    - Use gmsh.model.mesh.field.setNumber(mesh_field_id, "InField",
    ref_mesh_field_id): The "InField" option specifies that this new field will
    use the values from another existing field to compute its own values. In
    this case, it uses the same field with an identifier ref_mesh_field_id as a
    basis for computation. This means that the new field will be based on the
    values of another field with an identifier one less than the current 
    mesh_field_id. The "InField" option is useful for creating hierarchical mesh
    size fields, where one field depends on another, allowing more control over
    mesh refinement strategies.
    - Set the size constraints for the elements using 
    gmsh.model.mesh.field.setNumber(field_id, "SizeMin", size_min) and 
    gmsh.model.mesh.field.setNumber(field_id, "SizeMax", size_max).
    - Set the distance constraints for the elements using 
    gmsh.model.mesh.field.setNumber(field_id, "DistMin", dist_min) and
    gmsh.model.mesh.field.setNumber(field_id, "DistMax", dist_max).

    Args:
        mesh_field_id (int): A unique identifier for the threshold-based
          mesh size field.
        ref_mesh_field_id (int): The identifier of the existing mesh
          size field to base the new field on.
        size_min (float): Minimum size constraint for the elements.
        size_max (float): Maximum size constraint for the elements.
        dist_min (float): Minimum distance constraint for the elements.
        dist_max (float): Maximum distance constraint for the elements.
    """

    # Add a threshold-based mesh size field using the given field ID
    gmsh.model.mesh.field.add("Threshold", mesh_field_id)

    # Set the reference field that the new field will use to compute its
    #  values
    gmsh.model.mesh.field.setNumber(mesh_field_id, "InField", ref_mesh_field_id)

    # Set the size constraints for the elements.
    gmsh.model.mesh.field.setNumber(mesh_field_id, "SizeMin", size_min)
    gmsh.model.mesh.field.setNumber(mesh_field_id, "SizeMax", size_max)

    # Set the distance constraints for the elements
    gmsh.model.mesh.field.setNumber(mesh_field_id, "DistMin", dist_min)
    gmsh.model.mesh.field.setNumber(mesh_field_id, "DistMax", dist_max)


class GmshMesh:
    """Gmsh mesh.

    Attributes:
        geometry: A geometry_utils.GeometricCase object.
        global_refinement_factor (float): The refinement factor for global
          refinement of the mesh. A higher value results in a finer mesh
          overall, increasing the number of elements in the entire mesh, and
          leading to a more detailed representation of the geometry. Use this
          factor when you want to globally refine the mesh uniformly, without 
          specific local focus.
        local_refinement_factor (float): The refinement factor for local
          refinement of the mesh. This factor controls the local refinement
          level of the mesh and is typically used for refining specific regions
          or features of the mesh. A higher value for this factor indicates a
          finer mesh in the regions of interest, providing more detailed
          resolution around certain features. Use this factor when you want to
          focus on refining specific areas while keeping the rest of the mesh
          less refined.
    """

    def __init__(self,
                 geometry: geometry_utils.GeometricCase,
                 global_refinement_factor: Optional[float] = 1.0,
                 local_refinement_factor: Optional[float] = 10.0) -> None:
        """Initializes a GmshMesh object."""
        self.geometry = geometry
        self.global_refinement_factor = global_refinement_factor
        self.local_refinement_factor = local_refinement_factor

    def _create_mesh_with_gmsh(self) -> None:
        """Creates the mesh with Gmsh.

        To generate the mesh using Gmsh, we utilize mesh size fields: the 
        "Distance" and "Threshold" fields. Use add_mesh_field_distance to add
        the distance mesh field and add_mesh_field_threshold for the threshold
        mesh field in Gmsh.

        To control the element size using distance and threshold mesh size
        fields in Gmsh, we need to follow these steps:

        1. Define the distance-based mesh size field

        2. Define the threshold-based mesh size field

        3. Combine the distance and threshold-based fields:
            - To apply both distance and threshold fields, you can use the field 
            IDs of the previously defined fields and combine them using 
            mathematical operations (e.g., "Min" or "Max").
            - Use gmsh.model.mesh.field.add("Min", combined_mesh_field_id) to
            add a combined field that uses the minimum element size from the 
            distance and threshold fields.
            - Use gmsh.model.mesh.field.setNumbers(combined_mesh_field_id, 
            "FieldsList", list_mesh_field_id) to set the list of mesh fields 
            that will be combined in the newly created field with the 
            identifier combined_mesh_field_id.

        4. Set the mesh field as the background mesh field:
            - Use gmsh.model.mesh.field.setAsBackgroundMesh(
            combined_mesh_field_id) to set the combined field as the
            background mesh field, which will be used during mesh generation.

        5. Generate the mesh:
            - Finally, use gmsh.model.mesh.generate(dim=2) to generate the 2D
            mesh based on the defined mesh size fields and other meshing
            options.

        We can define both the global mesh size and the local mesh size (around
        the holes).
        """

        # Initialize the Gmsh API
        gmsh.initialize()

        # Disable Gmsh terminal output
        gmsh.option.setNumber("General.Terminal", 0)

        # Add a new model and set it as the current model
        gmsh.model.add("model")

        # Generate the plate with holes and get the boundary IDs
        (plate_curves_id, holes_curve_id
        ) = self.geometry.palte_with_holes_to_occ_and_get_boundary_ids()

        # Get mesh parameters
        (plate_mesh_offset, plate_predefined_element_size, holes_mesh_offset,
         holes_predefined_element_size) = self.geometry.get_mesh_params()

        # Add physical markers for the plate and holes
        gmsh.model.addPhysicalGroup(2, [1])

        # Calculate the maximum global mesh size
        global_mesh_size_max = (plate_predefined_element_size /
                                self.global_refinement_factor)

        mesh_field_id = 0
        mesh_field_threshol_id = []

        # Loop through plate boundaries to apply distance and threshold-based
        # mesh fields
        for id_plate_boundary in range(len(plate_curves_id)):

            mesh_field_id += 1
            add_mesh_field_distance(mesh_field_id,
                                    [plate_curves_id[id_plate_boundary]], 100)

            mesh_field_id += 1
            add_mesh_field_threshold(mesh_field_id, mesh_field_id - 1,
                                     global_mesh_size_max, global_mesh_size_max,
                                     0, plate_mesh_offset)
            mesh_field_threshol_id.append(mesh_field_id)

        # Check if local refinement is required
        if self.local_refinement_factor > 0.0:

            # Loop through holes to apply distance and threshold-based mesh
            # fields
            for id_hole in range(len(holes_curve_id)):

                # Calculate the local mesh size
                local_mesh_size = holes_predefined_element_size[
                    id_hole] / self.local_refinement_factor

                # Loop throught curves for each hole
                for hole_curve_id in holes_curve_id[id_hole]:

                    # Set the number of points to use along the hole's curve
                    # (100 for lines, 400 for other types)
                    if gmsh.model.getType(1, hole_curve_id) == "Line":
                        number_points = 100
                    else:
                        number_points = 400

                    mesh_field_id += 1
                    add_mesh_field_distance(mesh_field_id, [hole_curve_id],
                                            number_points)

                    mesh_field_id += 1
                    add_mesh_field_threshold(mesh_field_id, mesh_field_id - 1,
                                             local_mesh_size,
                                             global_mesh_size_max, 0,
                                             holes_mesh_offset[id_hole])
                    mesh_field_threshol_id.append(mesh_field_id)

        # Combine the distance and threshold-based fields
        combined_mesh_field_id = mesh_field_id + 1
        gmsh.model.mesh.field.add("Min", combined_mesh_field_id)
        gmsh.model.mesh.field.setNumbers(combined_mesh_field_id, "FieldsList",
                                         mesh_field_threshol_id)

        # Set the mesh field as the background mesh field
        gmsh.model.mesh.field.setAsBackgroundMesh(combined_mesh_field_id)

        # Set the meshing algorithm: Quasi-structured Quad method
        gmsh.option.setNumber("Mesh.Algorithm", 11)

        # Set mesh smoothing to improve mesh quality
        gmsh.option.set_number("Mesh.Smoothing", 1000)

        # Generate the mesh in 2D
        gmsh.model.mesh.generate(dim=2)

    def write_to_msh(self, msh_path: str) -> None:
        """Writes the GMSH mesh to MSH file.

        Args:
            msh_path (str): The mesh file path in MSH format.
        """
        self._create_mesh_with_gmsh()

        gmsh.write(msh_path)
        gmsh.finalize()
