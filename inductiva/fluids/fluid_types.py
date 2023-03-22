"""Physical properties of different fluid types."""

from dataclasses import dataclass


@dataclass(eq=True, frozen=True)
class FluidType:
    density: float
    kinematic_viscosity: float


WATER = FluidType(
    density=1e3,
    kinematic_viscosity=1e-6,
)

HONEY = FluidType(
    density=1360,
    kinematic_viscosity=7.36e-5,
)

OLIVE_OIL = FluidType(
    density=905,
    kinematic_viscosity=4.32e-5,
)

LIQUID_PROPANE = FluidType(
    density=580,
    kinematic_viscosity=2.73e-7,
)

JET_FUEL = FluidType(
    density=802.5,
    kinematic_viscosity=7.6e-6,
)

BEER = FluidType(
    density=1.007e3,
    kinematic_viscosity=1.8e-6,
)

GEAR_OIL = FluidType(
    density=860,
    kinematic_viscosity=4.2e-6,
)


def get_fluid_color(fluid):
    colors = {
        WATER: "blue",
        HONEY: "darkgoldenrod",
        OLIVE_OIL: "olivedrab",
        BEER: "goldenrod",
        GEAR_OIL: "peru",
        JET_FUEL: "darkgoldenrod",
        LIQUID_PROPANE: "blueviolet"
    }

    return colors[fluid]
