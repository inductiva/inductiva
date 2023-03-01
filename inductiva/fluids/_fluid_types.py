"""Physical properties of different fluid types."""

from inductiva_sph.sph_core.fluids import FluidProperties

WATER = FluidProperties(
    density = 1e3,
    kinematic_viscosity = 1e-6,
)

OLIVE_OIL = FluidProperties(
    density = 905,
    kinematic_viscosity = 4.3200e-5,
)

LIQUID_PROPANE = FluidProperties(
    density = 580,
    kinematic_viscosity = 2.73e-7,
)

JET_FUEL = FluidProperties(
    density = 802.5,
    kinematic_viscosity = 7.6e-6,
)

BEER = FluidProperties(
    density = 1.007e3,
    kinematic_viscosity = 1.8e-6,
)

GEAR_OIL = FluidProperties(
    density = 860,
    kinematic_viscosity = 4.2e-6,
)