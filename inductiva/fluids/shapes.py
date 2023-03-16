"""Classes defining arbitrary geometrical shapes."""

from dataclasses import dataclass, field
from typing import List


# Shapes.
@dataclass
class BaseShape:
    """Base shape."""
    pass


@dataclass
class Rectangle(BaseShape):
    """Rectangular shape."""
    dimensions: List[float] = field(default_factory=lambda: [1, 1])
    position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class Circle(BaseShape):
    """Circular shape."""
    radius: float = 1
    position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class Cube(BaseShape):
    """Cubic shape."""
    dimensions: List[float] = field(default_factory=lambda: [1, 1, 1])
    position: List[float] = field(default_factory=lambda: [0, 0, 0])


@dataclass
class Cylinder(BaseShape):
    """Cylindrical shape."""
    radius: float = 0.5
    height: float = 1
    position: List[float] = field(default_factory=lambda: [0, 0, 0])
