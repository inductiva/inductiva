"""Utils to create plates."""

from abc import ABC, abstractmethod
from typing import Tuple

import gmsh


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
    def get_plate_mesh_params(self):
        """Abstract method to get plate mesh parameters."""
        pass


class RectangularPlate(Plate):
    """Rectangular plate.

    Attributes:
        plate_type (str): Plate type.
        width (float): Plate width.
        length (float): Plate length.
    """

    def __init__(self, width: float, length: float) -> None:
        """Initializes a RectangularPlate object."""
        self.plate_type = "rectangular"
        self.width = width
        self.length = length

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
            "plate_type": "rectangular",
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

    def get_plate_mesh_params(self) -> Tuple[float, float]:
        """Gets the mesh parameters for the plate.

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
