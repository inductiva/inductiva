"""Classes that define a fluid tank scenario and simulate it via API."""

from dataclasses import dataclass
from typing import Optional

from inductiva_sph import sph_core

from inductiva.fluids.shapes import BaseShape
from inductiva.fluids.shapes import Rectangle
from inductiva.fluids.shapes import Circle
from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder
from inductiva.fluids._fluid_types import WATER


# Tank inlets.
@dataclass
class BaseTankInlet:
    """Base tank inlet."""
    fluid_velocity: float = 1


@dataclass
class RectangularTankInlet(BaseTankInlet, Rectangle):
    """Rectangular tank inlet."""
    pass


@dataclass
class CircularTankInlet(BaseTankInlet, Circle):
    """Circular tank inlet."""
    pass


# Tank outlets.
@dataclass
class BaseTankOutlet:
    """Base tank outlet."""
    pass


@dataclass
class CubicTankOutlet(BaseTankOutlet, Cube):
    """Cubic tank outlet."""
    pass


@dataclass
class CylindricalTankOutlet(BaseTankOutlet, Cylinder):
    """Cylindrical tank outlet."""
    pass


@dataclass
class FluidTank:
    """Fluid tank."""
    shape: BaseShape = Cube(position=[-0.5, -0.5, -0.5], dimensions=[1, 1, 1])
    fluid: sph_core.fluids.FluidProperties = WATER
    fluid_level: float = 0
    inlet: Optional[BaseTankInlet] = \
        CircularTankInlet(radius=0.1, position=[0, 0])
    outlet: Optional[BaseTankOutlet] = \
        CylindricalTankOutlet(radius=0.1, height=0.1, position=[0, 0, -0.1])

    def simulate(self):
        pass
