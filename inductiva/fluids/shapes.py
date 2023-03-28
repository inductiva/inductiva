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
