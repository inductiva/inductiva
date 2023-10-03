"""Utils to create plates."""

from abc import ABC, abstractmethod


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
