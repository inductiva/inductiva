"""Utils to create the geometric case."""

from typing import List, Optional, Tuple

import gmsh
import json

from inductiva.structures import holes, plates


def get_boundary_ids(entity_gmsh: int) -> List[int]:
    """Gets the IDs of the plate or hole boundaries.

    Args:
        entity_gmsh (int): The Gmsh entity ID representing either a plate or a
          hole's OpenCASCADE CAD representation.

    Returns:
        List[int]: A list of boundary IDs corresponding to the entity's
          boundaries.
    """
    boundaries = gmsh.model.getBoundary([(2, entity_gmsh)])

    return [[boundary_id[1]][0] for boundary_id in boundaries]


class GeometricCase:
    """Geometric case.

    The geometric case is characterized by a plate and a set of holes.

    Attributes:
        plate (RectangularPlate): Rectangular plate object.
        holes_list (List[Hole]): The holes objects.
    """

    def __init__(self,
                 plate: plates.RectangularPlate,
                 holes_list: Optional[List[holes.Hole]] = None) -> None:
        """Initializes a GeometricCase object."""
        self.plate = plate
        self.holes_list = holes_list

    def write_to_json(self, json_path: str) -> None:
        """Writes the geometric case to JSON file.

        Args:
            json_path (str): The JSON file path.
        """

        # Gemeometric case dictionary
        plate_dict = {"plate": self.plate.to_dict()}
        holes_dict = {"holes": [hole.to_dict() for hole in self.holes_list]}
        geom_case_dictionary = {**plate_dict, **holes_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(geom_case_dictionary, write_file, indent=4)

    def _get_holes_mesh_params(self) -> Tuple[List[float], List[float]]:
        """Gets the mesh generation parameters for all the hole.

        Metrics:
          - mesh_offset (float): Represents an offset for the curves of the
            holes, defining a region around the boundaries. Within this region,
            we have the ability to control the mesh elements size.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Returns:
            Tuple[List[float], List[float]]: A tuple containing two lists:
            - List of mesh offsets for each hole.
            - List of predefined element mesh sizes for each hole.
        """
        holes_mesh_offset = []
        holes_predefined_element_size = []

        for hole in self.holes_list:
            mesh_offset, predefined_element_size = hole.get_hole_mesh_params()

            holes_mesh_offset.append(mesh_offset)
            holes_predefined_element_size.append(predefined_element_size)

        return holes_mesh_offset, holes_predefined_element_size

    def _holes_to_occ_and_get_boundary_ids(
            self) -> Tuple[List[int], List[List[int]]]:
        """Converts list of holes to OpenCASCADE CAD, gets boundary IDs.

        Returns:
            Tuple[List[int], List[List[int]]]: A tuple containing two lists:
            - List of the Gmsh entity ID representing the hole's OpenCASCADE
            CAD representation for each hole.
            - List of lists, where each inner list contains the boundary IDs
            corresponding to the hole's boundaries for each hole.
        """

        holes_gmsh = []
        holes_boundary_ids = []

        for hole in self.holes_list:
            hole_gmsh = hole.to_occ()
            gmsh.model.occ.synchronize()

            hole_boundary_ids = get_boundary_ids(hole_gmsh)

            holes_gmsh.append(hole_gmsh)
            holes_boundary_ids.append(hole_boundary_ids)

        return holes_gmsh, holes_boundary_ids

    def palte_with_holes_to_occ_and_get_boundary_ids(
            self) -> Tuple[List[int], List[List[int]]]:
        """Converts palte with holes to OpenCASCADE CAD, gets boundary IDs.

        The process of generating the plate with holes is divided into 3 steps:

        1. Converts plate object to OpenCASCADE CAD representation and gets the
        IDs of the boundaries
        2. Converts holes objects to OpenCASCADE CAD representation and gets the
        IDs of the boundaries
        3. Removes holes from the plate in the OpenCASCADE CAD representation

        To remove the holes from the plate, the gmsh.model.occ.cut() function
        will be used.

        Returns:
            Tuple[List[int], List[List[int]]]:
                A tuple containing the following:
                - List of IDs of the plate's boundaries.
                - List of lists, where each inner list contains the IDs of the
                boundaries corresponding to the boundaries of each hole.
        """

        # 1. Converts the plate object to OpenCASCADE CAD representation and
        # gets the IDs of the boundaries
        plate_gmsh = self.plate.to_occ()
        gmsh.model.occ.synchronize()
        plate_boundary_ids = get_boundary_ids(plate_gmsh)

        # 2. Converts the hole objects to OpenCASCADE CAD representation and
        # gets the IDs of the boundaries
        (holes_gmsh,
         holes_boundary_ids) = self._holes_to_occ_and_get_boundary_ids()

        # 3. Removes holes from the plate in the OpenCASCADE CAD representation
        for hole_gmsh in holes_gmsh:
            gmsh.model.occ.cut([(2, plate_gmsh)], [(2, hole_gmsh)])
        gmsh.model.occ.synchronize()

        return plate_boundary_ids, holes_boundary_ids

    def get_mesh_params(self) -> Tuple[float, float, List[float], List[float]]:
        """Gets the mesh generation parameters for the palte with holes.

        Returns:
            Tuple[float, float, List[float], List[float]]:
                A tuple containing the following:
                - Mesh offset for the plate.
                - Predefined element mesh size for the plate.
                - List of mesh offsets for the holes.
                - List of predefined element sizes for the holes.
        """

        (plate_mesh_offset,
         plate_predefined_element_size) = self.plate.get_plate_mesh_params()
        (holes_mesh_offset,
         holes_predefined_element_size) = self._get_holes_mesh_params()

        return (plate_mesh_offset, plate_predefined_element_size,
                holes_mesh_offset, holes_predefined_element_size)
