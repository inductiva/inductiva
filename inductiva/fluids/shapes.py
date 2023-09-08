"""Classes defining arbitrary geometrical shapes."""

from abc import abstractmethod
from dataclasses import dataclass, field
from typing import List, Tuple


# Shapes.
@dataclass
class BaseShape:
    """Base shape."""
    pass

    @abstractmethod
    def get_bounding_box(self) -> Tuple[List[float], List[float]]:
        """Returns the bounding box of the shape."""
        pass

    @abstractmethod
    def to_dict(self) -> dict:
        """Returns a dictionary representation of the shape."""
        pass


@dataclass
class Base2DShape(BaseShape):
    """Base 2D shape."""
    position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class Base3DShape(BaseShape):
    """Base 3D shape."""
    position: List[float] = field(default_factory=lambda: [0, 0, 0])


@dataclass
class Rectangle(Base2DShape):
    """Rectangular shape."""
    dimensions: List[float] = field(default_factory=lambda: [1, 1])

    def get_bounding_box(self) -> Tuple[List[float], List[float]]:
        """Returns the bounding box of the rectangle."""
        return (
            self.position,
            [self.position[i] + self.dimensions[i] for i in range(2)],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the rectangle."""
        return {
            "type": "rectangle",
            "position": self.position,
            "dimensions": self.dimensions,
        }


@dataclass
class Circle(Base2DShape):
    """Circular shape."""
    radius: float = 1

    def get_bounding_box(self) -> Tuple[List[float], List[float]]:
        """Returns the bounding box of the circle."""
        return (
            [
                self.position[0] - self.radius,
                self.position[1] - self.radius,
            ],
            [
                self.position[0] + self.radius,
                self.position[1] + self.radius,
            ],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the circle."""
        return {
            "type": "circle",
            "position": self.position,
            "radius": self.radius,
        }


@dataclass
class Cube(Base3DShape):
    """Cubic shape."""
    dimensions: List[float] = field(default_factory=lambda: [1, 1, 1])

    def get_bounding_box(self) -> Tuple[List[float], List[float]]:
        """Returns the bounding box of the cube."""
        return (
            self.position,
            [self.position[i] + self.dimensions[i] for i in range(3)],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the cube."""
        return {
            "type": "cube",
            "position": self.position,
            "dimensions": self.dimensions,
        }


@dataclass
class Cylinder(Base3DShape):
    """Cylindrical shape."""
    radius: float = 0.5
    height: float = 1

    def get_bounding_box(self) -> Tuple[List[float], List[float]]:
        """Returns the bounding box of the cylinder."""
        return (
            [
                self.position[0] - self.radius,
                self.position[1] - self.radius,
                self.position[2],
            ],
            [
                self.position[0] + self.radius,
                self.position[1] + self.radius,
                self.position[2] + self.height,
            ],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the cylinder."""
        return {
            "type": "cylinder",
            "position": self.position,
            "radius": self.radius,
            "height": self.height,
        }
