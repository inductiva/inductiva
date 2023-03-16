"""Classes that define a fluid tank scenario and simulate it via API."""

from dataclasses import dataclass, field
from typing import List, Optional, Type

from inductiva_sph import sph_core

from inductiva.fluids._fluid_types import WATER


# Tank shapes.
@dataclass
class BaseTankShape:
    """Base tank shape."""
    pass


@dataclass
class CubicShape(BaseTankShape):
    """Cubic tank shape."""
    dimensions: List[float] = field(default_factory=lambda: [1, 1, 1])
    position: List[float] = field(default_factory=lambda: [-0.5, -0.5, -0.5])


@dataclass
class CylindricalShape(BaseTankShape):
    """Cylindrical tank shape."""
    radius: float = 0.5
    height: float = 1
    position: List[float] = field(default_factory=lambda: [0, 0, 0])


# Tank inlets.
@dataclass
class BaseTankInlet:
    """Base tank inlet."""
    fluid_velocity: float = 1


@dataclass
class CubicTankInlet(BaseTankInlet):
    """Cubic tank inlet."""
    dimensions: List[float] = field(default_factory=lambda: [0.1, 0.1])
    position: List[float] = field(default_factory=lambda: [-0.05, -0.05])


@dataclass
class CylindricalTankInlet(BaseTankInlet):
    """Cylindrical tank inlet."""
    radius: float = 0.1
    position: List[float] = field(default_factory=lambda: [0, 0])


# Tank outlets.
@dataclass
class BaseTankOutlet:
    """Base tank outlet."""
    pass


@dataclass
class CubicTankOutlet(BaseTankOutlet):
    """Cubic tank outlet."""
    dimensions: List[float] = field(default_factory=lambda: [0.1, 0.1, 0.1])
    position: List[float] = field(default_factory=lambda: [-0.05, -0.05])


@dataclass
class CylindricalTankOutlet(BaseTankOutlet):
    """Cylindrical tank outlet."""
    radius: float = 0.1
    height: float = 0.1
    position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class FluidTank:
    """Fluid tank."""
    shape: Type[BaseTankShape] = CubicShape()
    fluid: sph_core.fluids.FluidProperties = WATER
    fluid_level: float = 0
    inlet: Optional[Type[BaseTankInlet]] = CylindricalTankInlet()
    outlet: Optional[Type[BaseTankOutlet]] = CylindricalTankOutlet()
