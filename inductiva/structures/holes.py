"""Utils to create the holes."""

from abc import ABC, abstractmethod

import math


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


class CircularHole(Hole):
    """Circular hole.

    Attributes:
        radius (float): Hole radius.
    """

    def __init__(self, center_x: float, center_y: float, radius: float) -> None:
        """Initializes a CircularHole object."""
        super().__init__(center_x, center_y)
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
            "hole_type": "circular",
            "center_x": self.center_x,
            "center_y": self.center_y,
            "radius": self.radius
        }


class RectangularHole(Hole):
    """Rectangular hole.

    Attributes:
        half_size_x (float): Half size of the hole in the x-direction.
        half_size_y (float): Half size of the hole in the y-direction.
        angle (float): Positive angle of rotation in degrees around the hole
          center.
    """

    def __init__(self, center_x: float, center_y: float, half_size_x: float,
                 half_size_y: float, angle: float) -> None:
        """Initializes a RectangularHole object."""
        super().__init__(center_x, center_y)
        self.half_size_x = half_size_x
        self.half_size_y = half_size_y
        self.angle = angle

    def perimeter(self) -> float:
        """Calculate the perimeter of the hole.

        Returns:
            float: The calculated perimeter.
        """
        return self.half_size_x * 4 + self.half_size_y * 4

    def to_dict(self) -> dict:
        """Convert hole properties to a dictionary.

        Returns:
            dict: Hole properties.
        """
        return {
            "hole_type": "rectangular",
            "center_x": self.center_x,
            "center_y": self.center_y,
            "half_size_x": self.half_size_x,
            "half_size_y": self.half_size_y,
            "angle": self.angle
        }


class EllipticalHole(Hole):
    """Elliptical hole.

    Attributes:
        semi_axis_x (float): The semi-axis along the x-direction.
        semi_axis_y (float): The semi-axis along the y-direction.
        angle (float): Positive angle of rotation in degrees around the hole
          center.
    """

    def __init__(self, center_x: float, center_y: float, semi_axis_x: float,
                 semi_axis_y: float, angle: float) -> None:
        """Initializes a EllipticalHole object."""
        super().__init__(center_x, center_y)
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
            "hole_type": "elliptical",
            "center_x": self.center_x,
            "center_y": self.center_y,
            "semi_axis_x": self.semi_axis_x,
            "semi_axis_y": self.semi_axis_y,
            "angle": self.angle
        }
