"""Utils to create the geometric case."""

from abc import ABC, abstractmethod
from typing import List, Optional, Tuple

import gmsh
import json
import math


class Plate(ABC):
    """Abstract base class for plates."""

    @abstractmethod
    def perimeter(self):
        """Abstract method to calculate the perimeter of the plate."""
        pass

    @abstractmethod
    def to_dict(self):
        """Abstract method to convert the plate properties to a dictionary."""
        pass

    @abstractmethod
    def to_occ(self):
        """Abstract method to convert the plate to OpenCASCADE CAD."""
        pass

    @abstractmethod
    def get_boundary_ids(self):
        """Abstract method to get the IDs of the plate boundaries."""
        pass

    @abstractmethod
    def get_mesh_metrics(self):
        """Abstract method to get the metrics for mesh generation in Gmsh."""
        pass


class RectangularPlate(Plate):
    """Rectangular plate.

    Attributes:
        plate_type (str): Plate type.
        length (float): Plate length.
        width (float): Plate width.
    """

    def __init__(self, length: float, width: float) -> None:
        """Initializes a RectangularPlate object."""
        self.plate_type = "rectangular"
        self.length = length
        self.width = width

    def perimeter(self) -> float:
        """Calculate the perimeter of the plate.

        Returns:
            float: The calculated perimeter.
        """
        return 2 * (self.width + self.length)

    def to_dict(self) -> dict:
        """Convert plate properties to a dictionary.

        Returns:
            dict: Plate properties.
        """
        return {
            "plate_type": self.plate_type,
            "length": self.length,
            "width": self.width
        }

    def to_occ(self) -> int:
        """Converts plate to OpenCASCADE CAD representation.

        Returns:
          plate_gmsh (int): The Gmsh entity ID representing the plate's 
            OpenCASCADE CAD representation.
        """
        plate_gmsh = gmsh.model.occ.addRectangle(x=0,
                                                 y=0,
                                                 z=0,
                                                 dx=self.width,
                                                 dy=self.length)

        return plate_gmsh

    def get_boundary_ids(self, plate_gmsh: int) -> List[int]:
        """Gets the IDs of the plate boundaries.

        Args:
            plate_gmsh (int): The Gmsh entity ID representing the plate's 
              OpenCASCADE CAD representation.

        Returns:
            List[int]: A list of boundary IDs corresponding to the plate's 
              boundaries.
        """
        plate_boundaries = gmsh.model.getBoundary([(2, plate_gmsh)])

        return [[boundary_id[1]][0] for boundary_id in plate_boundaries]

    def get_mesh_metrics(self) -> Tuple[float, float]:
        """Gets the metrics for mesh generation in Gmsh.

        Metrics:
          - mesh_offset (float): Represents an offset for all the boundaries of 
            the plate, defining a region around the boundaries. Within this 
            region, we have the ability to control the mesh elements size.
            The offset is equal to half of the minimum size of the plate.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Returns:
          Tuple[float, float]: A tuple containing the mesh offset and the 
            predefined element mesh size for the plate.
      """
        mesh_offset = min(self.width, self.length) / 2
        predefined_element_size = self.perimeter() / 4

        return mesh_offset, predefined_element_size


class Hole(ABC):
    """Abstract base class for holes.

    Attributes:
        center_x (float): x-coordinate of the center.
        center_y (float): y-coordinate of the center.
    """

    def __init__(self, center_x: float, center_y: float) -> None:
        """Initializes a Hole object."""
        self.center_x = center_x
        self.center_y = center_y

    @abstractmethod
    def perimeter(self):
        """Abstract method to calculate the perimeter of the hole."""
        pass

    @abstractmethod
    def to_dict(self):
        """Abstract method to convert the hole properties to a dictionary."""
        pass

    @abstractmethod
    def to_occ(self):
        """Abstract method to convert the hole to OpenCASCADE CAD."""
        pass

    @abstractmethod
    def get_boundary_ids(self):
        """Abstract method to get the IDs of the hole boundaries."""
        pass

    @abstractmethod
    def get_mesh_metrics(self):
        """Abstract method to get the metrics for mesh generation in Gmsh."""
        pass


class CircularHole(Hole):
    """Circular hole.

    Attributes:
        hole_type (str): Hole type.
        radius (float): Hole radius.
    """

    def __init__(self, center_x: float, center_y: float, radius: float) -> None:
        """Initializes a CircularHole object."""
        super().__init__(center_x, center_y)
        self.hole_type = "circular"
        self.radius = radius

    def perimeter(self) -> float:
        """Calculate the perimeter of the hole.

        Returns:
            float: The calculated perimeter.
        """
        return 2 * math.pi * self.radius

    def to_dict(self) -> dict:
        """Convert hole properties to a dictionary.

        Returns:
            dict: Hole properties.
        """
        return {
            "hole_type": self.hole_type,
            "center_x": self.center_x,
            "center_y": self.center_y,
            "radius": self.radius
        }

    def to_occ(self) -> int:
        """Converts hole to OpenCASCADE CAD representation.

        Returns:
          hole_gmsh (int): The Gmsh entity ID representing the hole's 
            OpenCASCADE CAD representation.
        """
        hole_gmsh = gmsh.model.occ.addDisk(xc=self.center_x,
                                           yc=self.center_y,
                                           zc=0,
                                           rx=self.radius,
                                           ry=self.radius)
        return hole_gmsh

    def get_boundary_ids(self, hole_gmsh: int) -> List[int]:
        """Gets the IDs of the hole boundaries.

        Args:
            hole_gmsh (int): The Gmsh entity ID representing the hole's 
              OpenCASCADE CAD representation.

        Returns:
            List[int]: A list of boundary IDs corresponding to the hole's 
              boundaries.
        """
        hole_boundaries = gmsh.model.getBoundary([(2, hole_gmsh)])

        return [[boundary_id[1]][0] for boundary_id in hole_boundaries]

    def get_mesh_metrics(self) -> Tuple[float, float]:
        """Gets the metrics for mesh generation in Gmsh.

        Metrics:
          - mesh_offset (float): Represents an offset for the boundaries of the
            holes, defining a region around the boundaries. Within this region,
            we have the ability to control the mesh elements size. 
            The offset is equal to the radius.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Returns:
            Tuple[float, float]: The mesh offset and the predefined element mesh
              size for the hole.
        """
        mesh_offset = self.radius
        predefined_element_size = self.perimeter() / 4

        return mesh_offset, predefined_element_size


class RectangularHole(Hole):
    """Rectangular hole.

    Attributes:
        hole_type (str): Hole type.
        size_x (float): Size of the hole in the x-direction.
        size_y (float): Size of the hole in the y-direction.
        angle (float): Positive angle of rotation in degrees around the hole
          center.
    """

    def __init__(self, center_x: float, center_y: float, size_x: float,
                 size_y: float, angle: float) -> None:
        """Initializes a RectangularHole object."""
        super().__init__(center_x, center_y)
        self.hole_type = "rectangular"
        self.size_x = size_x
        self.size_y = size_y
        self.angle = angle

    def perimeter(self) -> float:
        """Calculate the perimeter of the hole.

        Returns:
            float: The calculated perimeter.
        """
        return self.size_x * 4 + self.size_y * 4

    def to_dict(self) -> dict:
        """Convert hole properties to a dictionary.

        Returns:
            dict: Hole properties.
        """
        return {
            "hole_type": self.hole_type,
            "center_x": self.center_x,
            "center_y": self.center_y,
            "size_x": self.size_x,
            "size_y": self.size_y,
            "angle": self.angle
        }

    def to_occ(self) -> int:
        """Converts hole to OpenCASCADE CAD representation.

        Returns:
          hole_gmsh (int): The Gmsh entity ID representing the hole's 
            OpenCASCADE CAD representation.
        """
        hole_gmsh = gmsh.model.occ.addRectangle(x=self.center_x - self.size_x,
                                                y=self.center_y - self.size_y,
                                                z=0,
                                                dx=self.size_x * 2,
                                                dy=self.size_y * 2)
        gmsh.model.occ.rotate(dimTags=[(2, hole_gmsh)],
                              x=self.center_x,
                              y=self.center_y,
                              z=0,
                              ax=0,
                              ay=0,
                              az=1,
                              angle=math.radians(self.angle))
        return hole_gmsh

    def get_boundary_ids(self, hole_gmsh: int) -> List[int]:
        """Gets the IDs of the hole boundaries.

        Args:
            hole_gmsh (int): The Gmsh entity ID representing the hole's 
              OpenCASCADE CAD representation.

        Returns:
            List[int]: A list of boundary IDs corresponding to the hole's 
              boundaries.
        """
        hole_boundaries = gmsh.model.getBoundary([(2, hole_gmsh)])

        return [[boundary_id[1]][0] for boundary_id in hole_boundaries]

    def get_mesh_metrics(self) -> Tuple[float, float]:
        """Gets the metrics for mesh generation in Gmsh.

        Metrics:
          - mesh_offset (float): Represents an offset for the boundaries of the
            holes, defining a region around the boundaries. Within this region,
            we have the ability to control the mesh elements size.
            The offset is equal to half of the minimum size of the hole.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Returns:
            Tuple[float, float]: The mesh offset and the predefined element mesh
              size for the hole.
        """
        mesh_offset = min(self.size_x, self.size_y) / 2
        predefined_element_size = self.perimeter() / 4

        return mesh_offset, predefined_element_size


class EllipticalHole(Hole):
    """Elliptical hole.

    Attributes:
        hole_type (str): Hole type.
        semi_axis_x (float): The semi-axis along the x-direction.
        semi_axis_y (float): The semi-axis along the y-direction.
        angle (float): Positive angle of rotation in degrees around the hole
          center.
    """

    def __init__(self, center_x: float, center_y: float, semi_axis_x: float,
                 semi_axis_y: float, angle: float) -> None:
        """Initializes a EllipticalHole object."""
        super().__init__(center_x, center_y)
        self.hole_type = "elliptical"
        self.semi_axis_x = semi_axis_x
        self.semi_axis_y = semi_axis_y
        self.angle = angle

    def perimeter(self) -> float:
        """Calculate the perimeter of the hole.

        Returns:
            float: The calculated perimeter.
        """
        return math.pi * (3 * (self.semi_axis_x + self.semi_axis_y) - math.sqrt(
            (3 * self.semi_axis_x + self.semi_axis_y) *
            (self.semi_axis_x + 3 * self.semi_axis_y)))

    def to_dict(self) -> dict:
        """Convert hole properties to a dictionary.

        Returns:
            dict: Hole properties.
        """

        return {
            "hole_type": self.hole_type,
            "center_x": self.center_x,
            "center_y": self.center_y,
            "semi_axis_x": self.semi_axis_x,
            "semi_axis_y": self.semi_axis_y,
            "angle": self.angle
        }

    def to_occ(self) -> int:
        """Converts hole to OpenCASCADE CAD representation.

        Returns:
            hole_gmsh (int): The Gmsh entity ID representing the hole's 
              OpenCASCADE CAD representation.
        """
        if self.semi_axis_x > self.semi_axis_y:
            hole_gmsh = gmsh.model.occ.addDisk(xc=self.center_x,
                                               yc=self.center_y,
                                               zc=0,
                                               rx=self.semi_axis_x,
                                               ry=self.semi_axis_y)
            gmsh.model.occ.rotate(dimTags=[(2, hole_gmsh)],
                                  x=self.center_x,
                                  y=self.center_y,
                                  z=0,
                                  ax=0,
                                  ay=0,
                                  az=1,
                                  angle=math.radians(self.angle))

        else:
            hole_gmsh = gmsh.model.occ.addDisk(xc=self.center_x,
                                               yc=self.center_y,
                                               zc=0,
                                               rx=self.semi_axis_y,
                                               ry=self.semi_axis_x)
            gmsh.model.occ.rotate(dimTags=[(2, hole_gmsh)],
                                  x=self.center_x,
                                  y=self.center_y,
                                  z=0,
                                  ax=0,
                                  ay=0,
                                  az=1,
                                  angle=math.radians(self.angle + 90))

        return hole_gmsh

    def get_boundary_ids(self, hole_gmsh: int) -> List[int]:
        """Gets the IDs of the hole boundaries.

        Args:
            hole_gmsh (int): The Gmsh entity ID representing the hole's 
              OpenCASCADE CAD representation.

        Returns:
            List[int]: A list of boundary IDs corresponding to the hole's 
              boundaries.
        """
        hole_boundaries = gmsh.model.getBoundary([(2, hole_gmsh)])

        return [[boundary_id[1]][0] for boundary_id in hole_boundaries]

    def get_mesh_metrics(self) -> Tuple[float, float]:
        """Gets the metrics for mesh generation in Gmsh.

        Metrics:
          - mesh_offset (float): Represents an offset for the boundaries of the
            holes, defining a region around the boundaries. Within this region,
            we have the ability to control the mesh elements size.
            The offset is equal to the minimum value semi-axis.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Returns:
            Tuple[float, float]: The mesh offset and the predefined element mesh
              size for the hole.
        """
        mesh_offset = min(self.semi_axis_x, self.semi_axis_y)
        predefined_element_size = self.perimeter / 4

        return mesh_offset, predefined_element_size


class GeometricCase:
    """Geometric case.

    The geometric case is characterized by a plate and a set of holes.

    Attributes:
        plate (RectangularPlate): Rectangular plate object.
        holes (List[Hole]): The holes objects.
    """

    def __init__(self,
                 plate: RectangularPlate,
                 holes: Optional[List[Hole]] = None) -> None:
        """Initializes a GeometricCase object."""
        self.plate = plate
        self.holes = holes

    def write_to_json(self, json_path: str) -> None:
        """Writes the geometric case to JSON file.

        Args:
            json_path (str): The JSON file path.
        """

        # Gemeometric case dictionary
        plate_dict = {"plate": self.plate.to_dict()}
        holes_dict = {"holes": [hole.to_dict() for hole in self.holes]}
        geom_case_dictionary = {**plate_dict, **holes_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(geom_case_dictionary, write_file, indent=4)

    def _holes_to_occ_and_get_boundary_ids(
            self) -> Tuple[List[int], List[List[int]]]:
        """Converts a list of hole objects to OpenCASCADE CAD representation and
        gets the IDs of the holes boundaries.

        Returns:
            Tuple[List[int], List[List[int]]]: A tuple containing two lists:
            - List of the Gmsh entity ID representing the hole's OpenCASCADE 
            CAD representation for each hole.
            - List of lists, where each inner list contains the boundary IDs 
            corresponding to the hole's boundaries for each hole.
        """

        holes_gmsh = []
        holes_boundary_ids = []

        for hole in self.holes:
            hole_gmsh = hole.to_occ()
            gmsh.model.occ.synchronize()

            hole_boundary_ids = hole.get_boundary_ids(hole_gmsh)

            holes_gmsh.append(hole_gmsh)
            holes_boundary_ids.append(hole_boundary_ids)

        return holes_gmsh, holes_boundary_ids

    def _get_holes_mesh_metrics(self) -> Tuple[List[float], List[float]]:
        """Gets the holes metrics for mesh generation in gmsh.

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

        for hole in self.holes:
            mesh_offset, predefined_element_size = hole.get_mesh_metrics()

            holes_mesh_offset.append(mesh_offset)
            holes_predefined_element_size.append(predefined_element_size)

        return holes_mesh_offset, holes_predefined_element_size

    def generate_plate_with_holes_get_boundary_ids_and_mesh_metrics(
        self
    ) -> Tuple[List[int], float, float, List[List[int]], List[float],
               List[float]]:
        """Generate the plate with holes in the OpenCASCADE CAD representation
         and gest the boundary IDs, and mesh metrics for the plate and holes.

        The process of generating the plate with holes is divided into 4vsteps:

        1. Converts plate object to OpenCASCADE CAD representation and gets the 
        IDs of the boundaries
        2. Converts holes objects to OpenCASCADE CAD representation and gets the 
        IDs of the boundaries
        3. Removes holes from the plate in the OpenCASCADE CAD representation
        4. Gets the holes and plate metrics required for building the mesh

        To remove the holes from the plate, the gmsh.model.occ.cut() function
        will be used.

        Returns:
            Tuple[List[int], float, float, List[List[int]], List[float], 
            List[float]]: 
                A tuple containing the following:
                - List of IDs of the plate's boundaries.
                - Mesh offset for the plate.
                - Predefined element mesh size for the plate.
                - List of lists, where each inner list contains the IDs of the
                boundaries corresponding to the boundaries of each hole.
                - List of mesh offsets for the holes.
                - List of predefined element sizes for the holes.
        """

        # 1. Converts the plate object to OpenCASCADE CAD representation and
        # gets the IDs of the boundaries
        plate_gmsh = self.plate.to_occ()
        gmsh.model.occ.synchronize()
        plate_boundary_ids = self.plate.get_boundary_ids(plate_gmsh)

        # 2. Converts the hole objects to OpenCASCADE CAD representation and
        # gets the IDs of the boundaries
        (holes_gmsh,
         holes_boundary_ids) = self._holes_to_occ_and_get_boundary_ids()

        # 3. Removes holes from the plate in the OpenCASCADE CAD representation
        for hole_gmsh in holes_gmsh:
            gmsh.model.occ.cut([(2, plate_gmsh)], [(2, hole_gmsh)])
        gmsh.model.occ.synchronize()

        # 4. Gets the holes and plate metrics required for building the mesh
        (plate_mesh_offset,
         plate_predefined_element_size) = self.plate.get_mesh_metrics()
        (holes_mesh_offset,
         holes_predefined_element_size) = self._get_holes_mesh_metrics()

        return (plate_boundary_ids, plate_mesh_offset,
                plate_predefined_element_size, holes_boundary_ids,
                holes_mesh_offset, holes_predefined_element_size)
