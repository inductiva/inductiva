"""Classes that define a fluid tank scenario and simulate it via API."""

from dataclasses import dataclass, field
from typing import List, Optional

from inductiva_sph import sph_core

from inductiva.fluids.shapes import BaseShape
from inductiva.fluids.shapes import Rectangle
from inductiva.fluids.shapes import Circle
from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder
from inductiva.fluids.fluid_types import WATER


# Tank inlets.
@dataclass
class BaseTankInlet:
    """Base tank inlet."""
    fluid_velocity: float = 1
    position: List[float] = field(default_factory=lambda: [0, 0])


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
    top_base_position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class CubicTankOutlet(BaseTankOutlet, Cube):
    """Cubic tank outlet."""
    pass


@dataclass
class CylindricalTankOutlet(BaseTankOutlet, Cylinder):
    """Cylindrical tank outlet."""
    pass


class FluidTank:
    """Fluid tank."""

    def __init__(
        self,
        shape: BaseShape = Cube(dimensions=[1, 1, 1]),
        fluid: sph_core.fluids.FluidProperties = WATER,
        fluid_level: float = 0,
        inlet: Optional[BaseTankInlet] = CircularTankInlet(radius=0.1,
                                                           position=[0, 0]),
        outlet: Optional[BaseTankOutlet] = CylindricalTankOutlet(
            radius=0.1, height=0.1, top_base_position=[0, 0]),
    ):
        self.shape = shape
        self.fluid = fluid
        self.fluid_level = fluid_level
        self.inlet = inlet
        self.outlet = outlet

    def simulate(self):
        pass
